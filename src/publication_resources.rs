//! Publication-associated external dataset catalog and download preparation.
//!
//! This module keeps literature-linked raw-data accessions in one portable
//! catalog and prepares deterministic manifests/download scripts without
//! forcing large downloads during routine resource-status checks.

use reqwest::blocking::Client;
use serde::{Deserialize, Serialize};
use sha1::{Digest, Sha1};
use std::{
    collections::BTreeSet,
    fs,
    io::{Read, Write},
    path::{Path, PathBuf},
};

pub const DEFAULT_PUBLICATION_RESOURCE_CATALOG_PATH: &str = "assets/publication_resources.json";
pub const DEFAULT_PUBLICATION_RESOURCE_CACHE_DIR: &str = "data/publication_resources";
const PUBLICATION_DATASET_LIST_SCHEMA: &str = "gentle.publication_dataset_list.v1";
const PUBLICATION_DATASET_STATUS_SCHEMA: &str = "gentle.publication_dataset_status.v1";
const PUBLICATION_DATASET_PREPARE_SCHEMA: &str = "gentle.publication_dataset_prepare.v1";
const PUBLICATION_DATASET_DOWNLOAD_MANIFEST_SCHEMA: &str =
    "gentle.publication_dataset_download_manifest.v1";
const BUILTIN_PUBLICATION_RESOURCE_CATALOG_JSON: &str =
    include_str!("../assets/publication_resources.json");

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PublicationResourceCatalog {
    pub schema: String,
    pub publication: PublicationRecord,
    pub datasets: Vec<PublicationDatasetCatalogEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PublicationRecord {
    pub id: String,
    pub title: String,
    pub doi: String,
    pub article_url: String,
    pub code_url: Option<String>,
    pub authors: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PublicationDatasetCatalogEntry {
    pub dataset_id: String,
    pub accession: String,
    pub repository: String,
    pub assay_kind: String,
    pub title: String,
    pub description: String,
    pub landing_url: String,
    pub api_url: Option<String>,
    pub file_listing_url: Option<String>,
    pub notes: Vec<String>,
    pub files: Vec<PublicationDatasetFile>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PublicationDatasetFile {
    pub file_name: String,
    pub url: String,
    pub category: String,
    pub size_label: Option<String>,
    pub size_bytes: Option<u64>,
    pub checksum_md5: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PublicationDatasetListEntry {
    pub dataset_id: String,
    pub accession: String,
    pub repository: String,
    pub assay_kind: String,
    pub title: String,
    pub description: String,
    pub landing_url: String,
    pub api_url: Option<String>,
    pub file_listing_url: Option<String>,
    pub declared_file_count: usize,
    pub categories: Vec<String>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PublicationDatasetListReport {
    pub schema: String,
    pub catalog_path: String,
    pub publication: PublicationRecord,
    pub filter: Option<String>,
    pub dataset_count: usize,
    pub datasets: Vec<PublicationDatasetListEntry>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PublicationDatasetFileStatus {
    pub file_name: String,
    pub url: String,
    pub category: String,
    pub size_label: Option<String>,
    pub expected_size_bytes: Option<u64>,
    pub expected_checksum_md5: Option<String>,
    pub local_path: String,
    pub exists: bool,
    pub local_size_bytes: Option<u64>,
    pub size_matches: Option<bool>,
    pub checksum_sha1: Option<String>,
    pub downloaded: bool,
    pub skipped: bool,
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PublicationDatasetStatusReport {
    pub schema: String,
    pub catalog_path: String,
    pub dataset_id: String,
    pub accession: String,
    pub repository: String,
    pub assay_kind: String,
    pub title: String,
    pub landing_url: String,
    pub api_url: Option<String>,
    pub file_listing_url: Option<String>,
    pub cache_dir: String,
    pub install_dir: String,
    pub declared_file_count: usize,
    pub local_file_count: usize,
    pub missing_file_count: usize,
    pub size_mismatch_count: usize,
    pub downloadable: bool,
    pub ready: bool,
    pub files: Vec<PublicationDatasetFileStatus>,
    pub notes: Vec<String>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PublicationDatasetPrepareReport {
    pub schema: String,
    pub catalog_path: String,
    pub dataset_id: String,
    pub accession: String,
    pub repository: String,
    pub cache_dir: String,
    pub install_dir: String,
    pub manifest_path: String,
    pub tsv_path: String,
    pub download_script_path: String,
    pub download_files: bool,
    pub max_files: Option<usize>,
    pub category_filters: Vec<String>,
    pub declared_file_count: usize,
    pub eligible_file_count: usize,
    pub planned_file_count: usize,
    pub downloaded_file_count: usize,
    pub skipped_file_count: usize,
    pub failed_file_count: usize,
    pub files: Vec<PublicationDatasetFileStatus>,
    pub notes: Vec<String>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PublicationResourceCollectionStatus {
    pub resource_id: String,
    pub display_name: String,
    pub support_status: String,
    pub catalog_path: String,
    pub cache_dir: String,
    pub catalog_valid: bool,
    pub dataset_count: usize,
    pub declared_file_count: usize,
    pub downloadable_dataset_count: usize,
    pub cached_dataset_count: usize,
    pub error: Option<String>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
struct PublicationDatasetDownloadManifest {
    schema: String,
    publication: PublicationRecord,
    dataset: PublicationDatasetListEntry,
    files: Vec<PublicationDatasetFile>,
    notes: Vec<String>,
}

pub fn list_publication_datasets(
    filter: Option<&str>,
    catalog_path: Option<&str>,
) -> Result<PublicationDatasetListReport, String> {
    let (catalog, label) = load_publication_resource_catalog(catalog_path)?;
    let filter = filter
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(str::to_string);
    let datasets: Vec<PublicationDatasetListEntry> = catalog
        .datasets
        .iter()
        .filter(|entry| dataset_matches_filter(entry, filter.as_deref()))
        .map(list_entry_for_dataset)
        .collect();
    Ok(PublicationDatasetListReport {
        schema: PUBLICATION_DATASET_LIST_SCHEMA.to_string(),
        catalog_path: label,
        publication: catalog.publication,
        filter,
        dataset_count: datasets.len(),
        datasets,
    })
}

pub fn publication_dataset_status(
    dataset_id: &str,
    catalog_path: Option<&str>,
    cache_dir: Option<&str>,
) -> Result<PublicationDatasetStatusReport, String> {
    let (catalog, label) = load_publication_resource_catalog(catalog_path)?;
    let entry = resolve_dataset(&catalog, dataset_id)?;
    let effective_cache_dir = effective_cache_dir(cache_dir);
    let install_dir = effective_cache_dir.join(sanitize_for_path(&entry.dataset_id));
    let files = file_statuses_for_entry(entry, &install_dir);
    let local_file_count = files.iter().filter(|row| row.exists).count();
    let missing_file_count = files.len().saturating_sub(local_file_count);
    let size_mismatch_count = files
        .iter()
        .filter(|row| row.size_matches == Some(false))
        .count();
    let warnings = unresolved_dataset_warnings(entry);
    Ok(PublicationDatasetStatusReport {
        schema: PUBLICATION_DATASET_STATUS_SCHEMA.to_string(),
        catalog_path: label,
        dataset_id: entry.dataset_id.clone(),
        accession: entry.accession.clone(),
        repository: entry.repository.clone(),
        assay_kind: entry.assay_kind.clone(),
        title: entry.title.clone(),
        landing_url: entry.landing_url.clone(),
        api_url: entry.api_url.clone(),
        file_listing_url: entry.file_listing_url.clone(),
        cache_dir: effective_cache_dir.display().to_string(),
        install_dir: install_dir.display().to_string(),
        declared_file_count: entry.files.len(),
        local_file_count,
        missing_file_count,
        size_mismatch_count,
        downloadable: !entry.files.is_empty(),
        ready: !entry.files.is_empty() && missing_file_count == 0 && size_mismatch_count == 0,
        files,
        notes: entry.notes.clone(),
        warnings,
    })
}

pub fn prepare_publication_dataset(
    dataset_id: &str,
    catalog_path: Option<&str>,
    cache_dir: Option<&str>,
    download_files: bool,
    max_files: Option<usize>,
    category_filters: &[String],
) -> Result<PublicationDatasetPrepareReport, String> {
    let (catalog, label) = load_publication_resource_catalog(catalog_path)?;
    let entry = resolve_dataset(&catalog, dataset_id)?.clone();
    let effective_cache_dir = effective_cache_dir(cache_dir);
    let install_dir = effective_cache_dir.join(sanitize_for_path(&entry.dataset_id));
    fs::create_dir_all(&install_dir).map_err(|e| {
        format!(
            "Could not create publication-resource directory '{}': {e}",
            install_dir.display()
        )
    })?;

    let category_filters = normalize_category_filters(category_filters);
    let eligible_files: Vec<PublicationDatasetFile> = entry
        .files
        .iter()
        .filter(|file| publication_file_matches_categories(file, &category_filters))
        .cloned()
        .collect();
    let planned_files: Vec<PublicationDatasetFile> = eligible_files
        .iter()
        .take(max_files.unwrap_or(usize::MAX))
        .cloned()
        .collect();
    let manifest_path = install_dir.join("manifest.json");
    let tsv_path = install_dir.join("download_manifest.tsv");
    let download_script_path = install_dir.join("download.sh");
    write_download_manifest_json(&manifest_path, &catalog.publication, &entry, &planned_files)?;
    write_download_manifest_tsv(&tsv_path, &planned_files)?;
    write_download_script(&download_script_path, &planned_files)?;

    let mut files = file_statuses_for_files(&planned_files, &install_dir);
    if download_files {
        let client = Client::builder()
            .build()
            .map_err(|e| format!("Could not build HTTP client for publication downloads: {e}"))?;
        for row in &mut files {
            if row.exists && row.size_matches != Some(false) {
                row.skipped = true;
                continue;
            }
            match download_publication_file(&client, row) {
                Ok(()) => {}
                Err(error) => row.error = Some(error),
            }
        }
    }

    let downloaded_file_count = files.iter().filter(|row| row.downloaded).count();
    let skipped_file_count = files.iter().filter(|row| row.skipped).count();
    let failed_file_count = files.iter().filter(|row| row.error.is_some()).count();
    let mut warnings = unresolved_dataset_warnings(&entry);
    if !download_files && !planned_files.is_empty() {
        warnings.push(
            "No bytes were downloaded; rerun with --download-files or execute download.sh to fetch the declared files.".to_string(),
        );
    }
    if max_files.is_some() && planned_files.len() < eligible_files.len() {
        warnings.push(format!(
            "Only {} of {} eligible files were planned because --max-files was set.",
            planned_files.len(),
            eligible_files.len()
        ));
    }
    if !category_filters.is_empty() {
        warnings.push(format!(
            "Applied category filter(s): {}.",
            category_filters.join(", ")
        ));
        if planned_files.is_empty() && !entry.files.is_empty() {
            warnings
                .push("No declared files matched the requested category filter(s).".to_string());
        }
    }
    Ok(PublicationDatasetPrepareReport {
        schema: PUBLICATION_DATASET_PREPARE_SCHEMA.to_string(),
        catalog_path: label,
        dataset_id: entry.dataset_id,
        accession: entry.accession,
        repository: entry.repository,
        cache_dir: effective_cache_dir.display().to_string(),
        install_dir: install_dir.display().to_string(),
        manifest_path: manifest_path.display().to_string(),
        tsv_path: tsv_path.display().to_string(),
        download_script_path: download_script_path.display().to_string(),
        download_files,
        max_files,
        category_filters,
        declared_file_count: entry.files.len(),
        eligible_file_count: eligible_files.len(),
        planned_file_count: planned_files.len(),
        downloaded_file_count,
        skipped_file_count,
        failed_file_count,
        files,
        notes: entry.notes,
        warnings,
    })
}

pub fn publication_resource_collection_status() -> PublicationResourceCollectionStatus {
    match load_publication_resource_catalog(None) {
        Ok((catalog, label)) => {
            let declared_file_count: usize =
                catalog.datasets.iter().map(|entry| entry.files.len()).sum();
            let downloadable_dataset_count = catalog
                .datasets
                .iter()
                .filter(|entry| !entry.files.is_empty())
                .count();
            let cache_root = PathBuf::from(DEFAULT_PUBLICATION_RESOURCE_CACHE_DIR);
            let cached_dataset_count = catalog
                .datasets
                .iter()
                .filter(|entry| {
                    cache_root
                        .join(sanitize_for_path(&entry.dataset_id))
                        .join("manifest.json")
                        .is_file()
                })
                .count();
            PublicationResourceCollectionStatus {
                resource_id: "publication_datasets".to_string(),
                display_name: "Publication-associated datasets".to_string(),
                support_status: if cached_dataset_count > 0 {
                    "catalog_with_prepared_manifests".to_string()
                } else {
                    "catalog_only".to_string()
                },
                catalog_path: label,
                cache_dir: DEFAULT_PUBLICATION_RESOURCE_CACHE_DIR.to_string(),
                catalog_valid: true,
                dataset_count: catalog.datasets.len(),
                declared_file_count,
                downloadable_dataset_count,
                cached_dataset_count,
                error: None,
                notes: vec![
                    "Use `resources list-publication-datasets` to inspect paper-linked accessions.".to_string(),
                    "Use `resources prepare-publication-dataset DATASET_ID` to write a manifest and download.sh without fetching large files.".to_string(),
                    "Add `--download-files` only when the local machine should fetch the declared raw files.".to_string(),
                ],
            }
        }
        Err(error) => PublicationResourceCollectionStatus {
            resource_id: "publication_datasets".to_string(),
            display_name: "Publication-associated datasets".to_string(),
            support_status: "catalog_error".to_string(),
            catalog_path: DEFAULT_PUBLICATION_RESOURCE_CATALOG_PATH.to_string(),
            cache_dir: DEFAULT_PUBLICATION_RESOURCE_CACHE_DIR.to_string(),
            catalog_valid: false,
            dataset_count: 0,
            declared_file_count: 0,
            downloadable_dataset_count: 0,
            cached_dataset_count: 0,
            error: Some(error),
            notes: vec![],
        },
    }
}

fn load_publication_resource_catalog(
    catalog_path: Option<&str>,
) -> Result<(PublicationResourceCatalog, String), String> {
    let (text, label) = if let Some(path) = catalog_path
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        (
            fs::read_to_string(path).map_err(|e| {
                format!("Could not read publication-resource catalog '{path}': {e}")
            })?,
            path.to_string(),
        )
    } else {
        (
            BUILTIN_PUBLICATION_RESOURCE_CATALOG_JSON.to_string(),
            DEFAULT_PUBLICATION_RESOURCE_CATALOG_PATH.to_string(),
        )
    };
    let catalog: PublicationResourceCatalog = serde_json::from_str(&text)
        .map_err(|e| format!("Could not parse publication-resource catalog '{label}': {e}"))?;
    if catalog.schema != "gentle.publication_datasets.v1" {
        return Err(format!(
            "Unsupported publication-resource catalog schema '{}' in '{}'",
            catalog.schema, label
        ));
    }
    Ok((catalog, label))
}

fn resolve_dataset<'a>(
    catalog: &'a PublicationResourceCatalog,
    query: &str,
) -> Result<&'a PublicationDatasetCatalogEntry, String> {
    let needle = query.trim();
    if needle.is_empty() {
        return Err("Publication dataset id/accession must not be empty".to_string());
    }
    catalog
        .datasets
        .iter()
        .find(|entry| {
            entry.dataset_id.eq_ignore_ascii_case(needle)
                || entry.accession.eq_ignore_ascii_case(needle)
        })
        .ok_or_else(|| format!("Publication dataset '{query}' was not found"))
}

fn dataset_matches_filter(entry: &PublicationDatasetCatalogEntry, filter: Option<&str>) -> bool {
    let Some(needle) = filter.map(|value| value.to_ascii_lowercase()) else {
        return true;
    };
    [
        entry.dataset_id.as_str(),
        entry.accession.as_str(),
        entry.repository.as_str(),
        entry.assay_kind.as_str(),
        entry.title.as_str(),
        entry.description.as_str(),
        entry.landing_url.as_str(),
    ]
    .into_iter()
    .chain(entry.notes.iter().map(String::as_str))
    .any(|value| value.to_ascii_lowercase().contains(&needle))
}

fn normalize_category_filters(values: &[String]) -> Vec<String> {
    values
        .iter()
        .map(|value| value.trim().to_ascii_lowercase())
        .filter(|value| !value.is_empty())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn publication_file_matches_categories(
    file: &PublicationDatasetFile,
    category_filters: &[String],
) -> bool {
    if category_filters.is_empty() {
        return true;
    }
    let category = file.category.trim().to_ascii_lowercase();
    category_filters.iter().any(|wanted| wanted == &category)
}

fn list_entry_for_dataset(entry: &PublicationDatasetCatalogEntry) -> PublicationDatasetListEntry {
    let categories: BTreeSet<String> = entry
        .files
        .iter()
        .map(|file| file.category.clone())
        .filter(|value| !value.trim().is_empty())
        .collect();
    PublicationDatasetListEntry {
        dataset_id: entry.dataset_id.clone(),
        accession: entry.accession.clone(),
        repository: entry.repository.clone(),
        assay_kind: entry.assay_kind.clone(),
        title: entry.title.clone(),
        description: entry.description.clone(),
        landing_url: entry.landing_url.clone(),
        api_url: entry.api_url.clone(),
        file_listing_url: entry.file_listing_url.clone(),
        declared_file_count: entry.files.len(),
        categories: categories.into_iter().collect(),
        notes: entry.notes.clone(),
    }
}

fn effective_cache_dir(cache_dir: Option<&str>) -> PathBuf {
    cache_dir
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(DEFAULT_PUBLICATION_RESOURCE_CACHE_DIR))
}

fn sanitize_for_path(value: &str) -> String {
    let sanitized: String = value
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || matches!(ch, '-' | '_' | '.') {
                ch
            } else {
                '_'
            }
        })
        .collect();
    if sanitized.is_empty() {
        "dataset".to_string()
    } else {
        sanitized
    }
}

fn file_statuses_for_entry(
    entry: &PublicationDatasetCatalogEntry,
    install_dir: &Path,
) -> Vec<PublicationDatasetFileStatus> {
    file_statuses_for_files(&entry.files, install_dir)
}

fn file_statuses_for_files(
    files: &[PublicationDatasetFile],
    install_dir: &Path,
) -> Vec<PublicationDatasetFileStatus> {
    files
        .iter()
        .map(|file| {
            let destination = install_dir.join(sanitize_file_name(&file.file_name));
            let local_size_bytes = destination.metadata().ok().map(|meta| meta.len());
            let exists = local_size_bytes.is_some();
            PublicationDatasetFileStatus {
                file_name: file.file_name.clone(),
                url: file.url.clone(),
                category: file.category.clone(),
                size_label: file.size_label.clone(),
                expected_size_bytes: file.size_bytes,
                expected_checksum_md5: file.checksum_md5.clone(),
                local_path: destination.display().to_string(),
                exists,
                local_size_bytes,
                size_matches: match (file.size_bytes, local_size_bytes) {
                    (Some(expected), Some(actual)) => Some(expected == actual),
                    _ => None,
                },
                checksum_sha1: None,
                downloaded: false,
                skipped: false,
                error: None,
            }
        })
        .collect()
}

fn sanitize_file_name(value: &str) -> String {
    let file_name = Path::new(value)
        .file_name()
        .and_then(|part| part.to_str())
        .unwrap_or(value);
    sanitize_for_path(file_name)
}

fn unresolved_dataset_warnings(entry: &PublicationDatasetCatalogEntry) -> Vec<String> {
    if entry.files.is_empty() {
        vec![format!(
            "No concrete file URLs are declared for {}; inspect {} or refresh the catalog once the archive exposes file links.",
            entry.accession, entry.landing_url
        )]
    } else {
        vec![]
    }
}

fn write_download_manifest_json(
    path: &Path,
    publication: &PublicationRecord,
    entry: &PublicationDatasetCatalogEntry,
    files: &[PublicationDatasetFile],
) -> Result<(), String> {
    let manifest = PublicationDatasetDownloadManifest {
        schema: PUBLICATION_DATASET_DOWNLOAD_MANIFEST_SCHEMA.to_string(),
        publication: publication.clone(),
        dataset: list_entry_for_dataset(entry),
        files: files.to_vec(),
        notes: entry.notes.clone(),
    };
    let mut text = serde_json::to_string_pretty(&manifest)
        .map_err(|e| format!("Could not serialize publication-resource manifest: {e}"))?;
    text.push('\n');
    fs::write(path, text).map_err(|e| {
        format!(
            "Could not write publication-resource manifest '{}': {e}",
            path.display()
        )
    })
}

fn write_download_manifest_tsv(
    path: &Path,
    files: &[PublicationDatasetFile],
) -> Result<(), String> {
    let mut text = String::from("file_name\tcategory\tsize_label\tsize_bytes\tchecksum_md5\turl\n");
    for file in files {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            file.file_name,
            file.category,
            file.size_label.as_deref().unwrap_or(""),
            file.size_bytes
                .map(|value| value.to_string())
                .unwrap_or_default(),
            file.checksum_md5.as_deref().unwrap_or(""),
            file.url
        ));
    }
    fs::write(path, text).map_err(|e| {
        format!(
            "Could not write publication-resource TSV manifest '{}': {e}",
            path.display()
        )
    })
}

fn write_download_script(path: &Path, files: &[PublicationDatasetFile]) -> Result<(), String> {
    let mut text =
        String::from("#!/usr/bin/env bash\nset -euo pipefail\ncd \"$(dirname \"$0\")\"\n");
    if files.is_empty() {
        text.push_str("echo 'No concrete file URLs are declared for this dataset yet.'\n");
    } else {
        for file in files {
            text.push_str(&format!(
                "curl -L --fail --retry 3 -o {} {}\n",
                shell_quote(&sanitize_file_name(&file.file_name)),
                shell_quote(&file.url)
            ));
        }
    }
    fs::write(path, text).map_err(|e| {
        format!(
            "Could not write publication-resource download script '{}': {e}",
            path.display()
        )
    })?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut permissions = fs::metadata(path)
            .map_err(|e| {
                format!(
                    "Could not inspect download script '{}': {e}",
                    path.display()
                )
            })?
            .permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(path, permissions).map_err(|e| {
            format!(
                "Could not mark download script '{}' executable: {e}",
                path.display()
            )
        })?;
    }
    Ok(())
}

fn shell_quote(value: &str) -> String {
    format!("'{}'", value.replace('\'', "'\"'\"'"))
}

fn download_publication_file(
    client: &Client,
    row: &mut PublicationDatasetFileStatus,
) -> Result<(), String> {
    let destination = PathBuf::from(&row.local_path);
    let mut response = client
        .get(&row.url)
        .send()
        .map_err(|e| format!("Could not download '{}': {e}", row.url))?;
    if !response.status().is_success() {
        return Err(format!(
            "Could not download '{}': HTTP {}",
            row.url,
            response.status()
        ));
    }
    let mut output = fs::File::create(&destination).map_err(|e| {
        format!(
            "Could not create publication-resource file '{}': {e}",
            destination.display()
        )
    })?;
    let mut hasher = Sha1::new();
    let mut buffer = [0u8; 64 * 1024];
    let mut bytes = 0u64;
    loop {
        let read = response
            .read(&mut buffer)
            .map_err(|e| format!("Could not read download '{}': {e}", row.url))?;
        if read == 0 {
            break;
        }
        output.write_all(&buffer[..read]).map_err(|e| {
            format!(
                "Could not write publication-resource file '{}': {e}",
                destination.display()
            )
        })?;
        hasher.update(&buffer[..read]);
        bytes = bytes.saturating_add(read as u64);
    }
    output.flush().map_err(|e| {
        format!(
            "Could not flush publication-resource file '{}': {e}",
            destination.display()
        )
    })?;
    row.exists = true;
    row.local_size_bytes = Some(bytes);
    row.size_matches = row.expected_size_bytes.map(|expected| expected == bytes);
    row.checksum_sha1 = Some(format!("{:x}", hasher.finalize()));
    if let Some(expected) = row.expected_size_bytes
        && expected != bytes
    {
        return Err(format!(
            "Downloaded '{}' size was {bytes} bytes, expected {expected} bytes",
            row.file_name
        ));
    }
    row.downloaded = true;
    Ok(())
}
