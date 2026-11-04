//! Interpretation-local FASTA index reuse with fresh resolution and final checks.
//! No file handles, failed lookups, sequence windows or global state are cached.

use super::*;
use std::time::SystemTime;

/// One interpretation's index reuse; validate before returning a report.
/// Catalog resolution and manifest/file admission still run for every window.
pub(crate) struct PreparedSequenceRegionReader<'a> {
    catalog: &'a GenomeCatalog,
    cache: Option<SequenceIndexCache>,
}

impl<'a> PreparedSequenceRegionReader<'a> {
    pub(crate) fn new(catalog: &'a GenomeCatalog) -> Self {
        Self {
            catalog,
            cache: Some(SequenceIndexCache::default()),
        }
    }

    #[cfg(test)]
    pub(crate) fn uncached(catalog: &'a GenomeCatalog) -> Self {
        Self {
            catalog,
            cache: None,
        }
    }

    pub(crate) fn get_sequence_region_with_cache(
        &mut self,
        genome_id: &str,
        chromosome: &str,
        start_1based: usize,
        end_1based: usize,
        cache_dir_override: Option<&str>,
    ) -> Result<String, String> {
        self.catalog.read_prepared_sequence_region(
            genome_id,
            chromosome,
            start_1based,
            end_1based,
            cache_dir_override,
            true,
            self.cache.as_mut(),
        )
    }

    /// Rehash bound manifests/indexes, including same-size/restored-mtime edits.
    /// This supplements, never replaces, the caller's BLAST reference checks.
    pub(crate) fn validate_unchanged(&mut self) -> Result<(), String> {
        if let Some(cache) = &mut self.cache {
            cache.validate_unchanged()?;
        }
        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn index_load_count(&self) -> usize {
        self.cache.as_ref().map_or(0, |cache| cache.index_loads)
    }
}

#[derive(PartialEq, Eq)]
struct IndexStamp {
    length: u64,
    modified: Option<SystemTime>,
    created: Option<SystemTime>,
}

impl IndexStamp {
    fn read(path: &Path) -> Result<Self, String> {
        let metadata = fs::metadata(path)
            .map_err(|error| format!("Could not open FASTA index '{}': {error}", path.display()))?;
        Ok(Self {
            length: metadata.len(),
            modified: metadata.modified().ok(),
            created: metadata.created().ok(),
        })
    }
}

struct CachedIndex {
    entries: HashMap<String, FastaIndexEntry>,
    sha256: String,
    stamp: IndexStamp,
}

#[derive(Default)]
pub(super) struct SequenceIndexCache {
    indexes: BTreeMap<PathBuf, CachedIndex>,
    manifests: BTreeMap<PathBuf, (String, GenomeInstallManifest)>,
    invalid: Option<String>,
    #[cfg(test)]
    index_loads: usize,
}

fn index_fingerprint(path: &Path) -> Result<String, String> {
    // Keep the legacy open/read error wording for ordinary unavailable inputs.
    drop(
        File::open(path)
            .map_err(|error| format!("Could not open FASTA index '{}': {error}", path.display()))?,
    );
    crate::digest_utils::sha256_file_hex(path)
        .map_err(|error| format!("Could not read FASTA index '{}': {error}", path.display()))
}

fn manifest_fingerprint(path: &Path) -> Result<String, String> {
    crate::digest_utils::sha256_file_hex(path).map_err(|error| {
        format!(
            "Could not recheck prepared manifest '{}': {error}",
            path.display()
        )
    })
}

fn changed(path: &Path) -> String {
    format!(
        "Prepared reference changed during interpretation: '{}'",
        path.display()
    )
}

impl SequenceIndexCache {
    pub(super) fn bind_manifest(
        &mut self,
        path: &Path,
        manifest: &GenomeInstallManifest,
    ) -> Result<(), String> {
        if let Some(error) = &self.invalid {
            return Err(error.clone());
        }
        let digest = manifest_fingerprint(path)?;
        if let Some((expected, _)) = self.manifests.get(path) {
            if &digest != expected {
                let error = changed(path);
                self.invalid = Some(error.clone());
                return Err(error);
            }
        } else {
            self.manifests
                .insert(path.to_path_buf(), (digest, manifest.clone()));
        }
        Ok(())
    }

    pub(super) fn index(
        &mut self,
        path: &Path,
    ) -> Result<&HashMap<String, FastaIndexEntry>, String> {
        if let Some(error) = &self.invalid {
            return Err(error.clone());
        }
        if let Some(cached) = self.indexes.get_mut(path) {
            let check = (|| {
                let stamp = IndexStamp::read(path)?;
                if stamp != cached.stamp {
                    if index_fingerprint(path)? != cached.sha256 {
                        return Err(changed(path));
                    }
                    cached.stamp = stamp;
                }
                Ok(())
            })();
            if let Err(error) = check {
                self.invalid = Some(error.clone());
                return Err(error);
            }
        } else {
            let before = index_fingerprint(path)?;
            let entries = load_fasta_index(path)?;
            let stamp = IndexStamp::read(path)?;
            if index_fingerprint(path)? != before {
                let error = changed(path);
                self.invalid = Some(error.clone());
                return Err(error);
            }
            self.indexes.insert(
                path.to_path_buf(),
                CachedIndex {
                    entries,
                    sha256: before,
                    stamp,
                },
            );
            #[cfg(test)]
            {
                self.index_loads += 1;
            }
        }
        Ok(&self
            .indexes
            .get(path)
            .expect("successful index load")
            .entries)
    }

    fn validate_unchanged(&mut self) -> Result<(), String> {
        if let Some(error) = &self.invalid {
            return Err(error.clone());
        }
        let result = (|| {
            for (path, (expected, manifest)) in &self.manifests {
                GenomeCatalog::validate_manifest_files(manifest)?;
                if manifest_fingerprint(path)? != *expected {
                    return Err(changed(path));
                }
            }
            // Always rehash: equal size/timestamps are not evidence of equal bytes.
            for (path, cached) in &self.indexes {
                if index_fingerprint(path)? != cached.sha256 {
                    return Err(changed(path));
                }
            }
            Ok(())
        })();
        if let Err(error) = &result {
            self.invalid = Some(error.clone());
        }
        result
    }
}
