use crate::feature_location::feature_is_reverse;
use flate2::read::GzDecoder;
use reqwest::blocking::get;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

pub const DEFAULT_GENOME_CATALOG_PATH: &str = "assets/genomes.json";
pub const DEFAULT_HELPER_GENOME_CATALOG_PATH: &str = "assets/helper_genomes.json";
pub const DEFAULT_GENOME_CACHE_DIR: &str = "data/genomes";

/// Catalog entry describing where to fetch one genome assembly and annotation.
#[derive(Default, Deserialize, Serialize, Debug, Clone)]
pub struct GenomeCatalogEntry {
    pub ncbi_taxonomy_id: Option<u32>,
    pub description: Option<String>,
    pub ncbi_assembly_accession: Option<String>,
    pub ncbi_assembly_name: Option<String>,
    pub genbank_accession: Option<String>,
    pub genbank_reference: Option<String>,
    pub local_variant_unpublished: Option<bool>,
    pub sequence_remote: Option<String>,
    pub annotations_remote: Option<String>,
    pub sequence_local: Option<String>,
    pub annotations_local: Option<String>,
    #[serde(default = "default_cache_dir")]
    pub cache_dir: Option<String>,
}

fn default_cache_dir() -> Option<String> {
    Some(DEFAULT_GENOME_CACHE_DIR.to_string())
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct GenomeInstallManifest {
    genome_id: String,
    sequence_source: String,
    annotation_source: String,
    sequence_path: String,
    annotation_path: String,
    fasta_index_path: String,
    gene_index_path: Option<String>,
    installed_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrepareGenomeReport {
    pub genome_id: String,
    pub reused_existing: bool,
    pub sequence_path: String,
    pub annotation_path: String,
    pub fasta_index_path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrepareGenomeProgress {
    pub genome_id: String,
    pub phase: String,
    pub item: String,
    pub bytes_done: u64,
    pub bytes_total: Option<u64>,
    pub percent: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomeGeneRecord {
    pub chromosome: String,
    pub start_1based: usize,
    pub end_1based: usize,
    pub strand: Option<char>,
    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub biotype: Option<String>,
}

#[derive(Debug, Clone)]
struct FastaIndexEntry {
    length: u64,
    offset: u64,
    line_bases: u64,
    line_bytes: u64,
}

#[derive(Debug, Clone, Default)]
pub struct GenomeCatalog {
    entries: HashMap<String, GenomeCatalogEntry>,
    catalog_base_dir: PathBuf,
}

impl GenomeCatalog {
    pub fn from_json_file(path: &str) -> Result<Self, String> {
        let text = fs::read_to_string(path)
            .map_err(|e| format!("Could not read genome catalog '{path}': {e}"))?;
        let entries: HashMap<String, GenomeCatalogEntry> = serde_json::from_str(&text)
            .map_err(|e| format!("Could not parse genome catalog '{path}': {e}"))?;
        let base = Path::new(path)
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        Ok(Self {
            entries,
            catalog_base_dir: base,
        })
    }

    pub fn list_genomes(&self) -> Vec<String> {
        let mut names: Vec<String> = self.entries.keys().cloned().collect();
        names.sort_unstable();
        names
    }

    pub fn is_prepared(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<bool, String> {
        let entry = self.entry(genome_id)?;
        let install_dir = self.install_dir(genome_id, entry, cache_dir_override);
        let manifest_path = install_dir.join("manifest.json");
        if !manifest_path.exists() {
            return Ok(false);
        }
        let manifest = Self::load_manifest(&manifest_path)?;
        if Self::validate_manifest_files(&manifest).is_err() {
            return Ok(false);
        }
        let gene_index_path = manifest
            .gene_index_path
            .as_ref()
            .map(PathBuf::from)
            .unwrap_or_else(|| install_dir.join("genes.json"));
        Ok(gene_index_path.exists())
    }

    pub fn prepare_genome_once(&self, genome_id: &str) -> Result<PrepareGenomeReport, String> {
        self.prepare_genome_once_with_cache(genome_id, None)
    }

    pub fn prepare_genome_once_with_cache(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<PrepareGenomeReport, String> {
        let mut noop = |_| {};
        self.prepare_genome_once_with_progress(genome_id, cache_dir_override, &mut noop)
    }

    pub fn prepare_genome_once_with_progress(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress),
    ) -> Result<PrepareGenomeReport, String> {
        let entry = self.entry(genome_id)?;
        let install_dir = self.install_dir(genome_id, entry, cache_dir_override);
        fs::create_dir_all(&install_dir).map_err(|e| {
            format!(
                "Could not create genome cache dir '{}': {e}",
                install_dir.display()
            )
        })?;
        let manifest_path = install_dir.join("manifest.json");

        if manifest_path.exists() {
            let mut manifest = Self::load_manifest(&manifest_path)?;
            Self::validate_manifest_files(&manifest)?;
            let gene_index_path = manifest
                .gene_index_path
                .as_ref()
                .map(PathBuf::from)
                .unwrap_or_else(|| install_dir.join("genes.json"));
            if !gene_index_path.exists() {
                build_gene_index_file(
                    Path::new(&manifest.annotation_path),
                    &gene_index_path,
                    |done, total| {
                        on_progress(PrepareGenomeProgress {
                            genome_id: genome_id.to_string(),
                            phase: "index_genes".to_string(),
                            item: canonical_or_display(Path::new(&manifest.annotation_path)),
                            bytes_done: done,
                            bytes_total: total,
                            percent: total.and_then(|t| {
                                if t == 0 {
                                    None
                                } else {
                                    Some((done as f64 / t as f64) * 100.0)
                                }
                            }),
                        });
                    },
                )?;
                manifest.gene_index_path = Some(canonical_or_display(&gene_index_path));
                Self::write_manifest(&manifest_path, &manifest)?;
            }
            on_progress(PrepareGenomeProgress {
                genome_id: genome_id.to_string(),
                phase: "ready".to_string(),
                item: "manifest".to_string(),
                bytes_done: 0,
                bytes_total: None,
                percent: Some(100.0),
            });
            return Ok(PrepareGenomeReport {
                genome_id: genome_id.to_string(),
                reused_existing: true,
                sequence_path: manifest.sequence_path,
                annotation_path: manifest.annotation_path,
                fasta_index_path: manifest.fasta_index_path,
            });
        }

        let sequence_source = self.resolve_source(
            genome_id,
            "sequence",
            entry.sequence_local.as_ref(),
            entry.sequence_remote.as_ref(),
            entry,
        )?;
        let annotation_source = self.resolve_source(
            genome_id,
            "annotation",
            entry.annotations_local.as_ref(),
            entry.annotations_remote.as_ref(),
            entry,
        )?;

        let sequence_path = install_dir.join("sequence.fa");
        let annotation_ext = infer_annotation_extension(&annotation_source);
        let annotation_path = install_dir.join(format!("annotation.{annotation_ext}"));
        let fasta_index_path = install_dir.join("sequence.fa.fai");
        let gene_index_path = install_dir.join("genes.json");

        materialize_source_with_progress(&sequence_source, &sequence_path, |done, total| {
            on_progress(PrepareGenomeProgress {
                genome_id: genome_id.to_string(),
                phase: "download_sequence".to_string(),
                item: sequence_source.clone(),
                bytes_done: done,
                bytes_total: total,
                percent: total.and_then(|t| {
                    if t == 0 {
                        None
                    } else {
                        Some((done as f64 / t as f64) * 100.0)
                    }
                }),
            });
        })?;
        materialize_source_with_progress(&annotation_source, &annotation_path, |done, total| {
            on_progress(PrepareGenomeProgress {
                genome_id: genome_id.to_string(),
                phase: "download_annotation".to_string(),
                item: annotation_source.clone(),
                bytes_done: done,
                bytes_total: total,
                percent: total.and_then(|t| {
                    if t == 0 {
                        None
                    } else {
                        Some((done as f64 / t as f64) * 100.0)
                    }
                }),
            });
        })?;
        build_fasta_index_with_progress(&sequence_path, &fasta_index_path, |done, total| {
            on_progress(PrepareGenomeProgress {
                genome_id: genome_id.to_string(),
                phase: "index_fasta".to_string(),
                item: canonical_or_display(&sequence_path),
                bytes_done: done,
                bytes_total: total,
                percent: total.and_then(|t| {
                    if t == 0 {
                        None
                    } else {
                        Some((done as f64 / t as f64) * 100.0)
                    }
                }),
            });
        })?;
        build_gene_index_file(&annotation_path, &gene_index_path, |done, total| {
            on_progress(PrepareGenomeProgress {
                genome_id: genome_id.to_string(),
                phase: "index_genes".to_string(),
                item: canonical_or_display(&annotation_path),
                bytes_done: done,
                bytes_total: total,
                percent: total.and_then(|t| {
                    if t == 0 {
                        None
                    } else {
                        Some((done as f64 / t as f64) * 100.0)
                    }
                }),
            });
        })?;

        let manifest = GenomeInstallManifest {
            genome_id: genome_id.to_string(),
            sequence_source,
            annotation_source,
            sequence_path: canonical_or_display(&sequence_path),
            annotation_path: canonical_or_display(&annotation_path),
            fasta_index_path: canonical_or_display(&fasta_index_path),
            gene_index_path: Some(canonical_or_display(&gene_index_path)),
            installed_at_unix_ms: now_unix_ms(),
        };
        Self::write_manifest(&manifest_path, &manifest)?;

        on_progress(PrepareGenomeProgress {
            genome_id: genome_id.to_string(),
            phase: "ready".to_string(),
            item: "manifest".to_string(),
            bytes_done: 0,
            bytes_total: None,
            percent: Some(100.0),
        });

        Ok(PrepareGenomeReport {
            genome_id: genome_id.to_string(),
            reused_existing: false,
            sequence_path: manifest.sequence_path,
            annotation_path: manifest.annotation_path,
            fasta_index_path: manifest.fasta_index_path,
        })
    }

    pub fn get_sequence_region(
        &self,
        genome_id: &str,
        chromosome: &str,
        start_1based: usize,
        end_1based: usize,
    ) -> Result<String, String> {
        self.get_sequence_region_with_cache(genome_id, chromosome, start_1based, end_1based, None)
    }

    pub fn get_sequence_region_with_cache(
        &self,
        genome_id: &str,
        chromosome: &str,
        start_1based: usize,
        end_1based: usize,
        cache_dir_override: Option<&str>,
    ) -> Result<String, String> {
        if start_1based == 0 {
            return Err("Coordinates must be 1-based (start >= 1)".to_string());
        }
        if end_1based < start_1based {
            return Err(format!(
                "Invalid interval: start ({start_1based}) is greater than end ({end_1based})"
            ));
        }

        let entry = self.entry(genome_id)?;
        let manifest_path = self
            .install_dir(genome_id, entry, cache_dir_override)
            .join("manifest.json");
        if !manifest_path.exists() {
            return Err(format!(
                "Genome '{genome_id}' is not prepared locally. Run PrepareGenome first."
            ));
        }

        let manifest = Self::load_manifest(&manifest_path)?;
        Self::validate_manifest_files(&manifest)?;

        let index = load_fasta_index(Path::new(&manifest.fasta_index_path))?;
        let mut chr_candidates = Vec::new();
        chr_candidates.push(chromosome.to_string());
        if let Some(trimmed) = chromosome.strip_prefix("chr") {
            chr_candidates.push(trimmed.to_string());
        } else {
            chr_candidates.push(format!("chr{chromosome}"));
        }

        let index_entry = chr_candidates
            .iter()
            .find_map(|name| index.get(name))
            .ok_or_else(|| {
                format!(
                    "Chromosome/contig '{}' not found in genome '{}'",
                    chromosome, genome_id
                )
            })?;

        let end_u64 = end_1based as u64;
        if end_u64 > index_entry.length {
            return Err(format!(
                "Requested end {} exceeds chromosome length {}",
                end_1based, index_entry.length
            ));
        }

        read_fasta_slice(
            Path::new(&manifest.sequence_path),
            index_entry,
            start_1based as u64,
            end_1based as u64,
        )
    }

    pub fn list_gene_regions(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<Vec<GenomeGeneRecord>, String> {
        let entry = self.entry(genome_id)?;
        let install_dir = self.install_dir(genome_id, entry, cache_dir_override);
        let manifest_path = install_dir.join("manifest.json");
        if !manifest_path.exists() {
            return Err(format!(
                "Genome '{genome_id}' is not prepared locally. Run PrepareGenome first."
            ));
        }

        let mut manifest = Self::load_manifest(&manifest_path)?;
        Self::validate_manifest_files(&manifest)?;
        let gene_index_path = manifest
            .gene_index_path
            .as_ref()
            .map(PathBuf::from)
            .unwrap_or_else(|| install_dir.join("genes.json"));
        if !gene_index_path.exists() {
            build_gene_index_file(
                Path::new(&manifest.annotation_path),
                &gene_index_path,
                |_done, _total| {},
            )?;
            manifest.gene_index_path = Some(canonical_or_display(&gene_index_path));
            Self::write_manifest(&manifest_path, &manifest)?;
        }
        load_gene_index_file(&gene_index_path)
    }

    fn entry(&self, genome_id: &str) -> Result<&GenomeCatalogEntry, String> {
        self.entries
            .get(genome_id)
            .ok_or_else(|| format!("Genome '{genome_id}' is not present in the catalog"))
    }

    fn install_dir(
        &self,
        genome_id: &str,
        entry: &GenomeCatalogEntry,
        cache_dir_override: Option<&str>,
    ) -> PathBuf {
        let base = cache_dir_override
            .map(|raw| self.resolve_local_path(raw))
            .or_else(|| {
                entry
                    .cache_dir
                    .as_ref()
                    .map(|raw| self.resolve_local_path(raw))
            })
            .unwrap_or_else(|| PathBuf::from(DEFAULT_GENOME_CACHE_DIR));
        base.join(sanitize_for_path(genome_id))
    }

    fn resolve_source(
        &self,
        genome_id: &str,
        kind: &str,
        local: Option<&String>,
        remote: Option<&String>,
        entry: &GenomeCatalogEntry,
    ) -> Result<String, String> {
        let mut missing_local_path: Option<PathBuf> = None;
        if let Some(local_raw) = local {
            let local_path = self.resolve_local_path(local_raw);
            if local_path.exists() {
                return Ok(canonical_or_display(&local_path));
            }
            missing_local_path = Some(local_path);
        }
        if let Some(remote_src) = remote {
            return Ok(remote_src.clone());
        }
        if let Some(ncbi_source) = resolve_ncbi_assembly_source(kind, entry)? {
            return Ok(ncbi_source);
        }
        if let Some(genbank_source) = resolve_genbank_accession_source(kind, entry)? {
            return Ok(genbank_source);
        }
        if let Some(local_path) = missing_local_path {
            return Err(format!(
                "Genome '{genome_id}' has a {kind}_local path '{}', but that file does not exist and no remote/NCBI fallback is configured",
                local_path.display()
            ));
        }
        Err(format!(
            "Genome '{genome_id}' does not provide {kind}_local/{kind}_remote and has no resolvable NCBI assembly/GenBank accession source"
        ))
    }

    fn resolve_local_path(&self, raw: &str) -> PathBuf {
        if let Some(stripped) = raw.strip_prefix("file://") {
            return PathBuf::from(stripped);
        }
        let p = Path::new(raw);
        if p.is_absolute() {
            p.to_path_buf()
        } else {
            self.catalog_base_dir.join(p)
        }
    }

    fn load_manifest(path: &Path) -> Result<GenomeInstallManifest, String> {
        let text = fs::read_to_string(path)
            .map_err(|e| format!("Could not read genome manifest '{}': {e}", path.display()))?;
        serde_json::from_str(&text)
            .map_err(|e| format!("Could not parse genome manifest '{}': {e}", path.display()))
    }

    fn write_manifest(path: &Path, manifest: &GenomeInstallManifest) -> Result<(), String> {
        let text = serde_json::to_string_pretty(manifest)
            .map_err(|e| format!("Could not serialize genome manifest: {e}"))?;
        fs::write(path, text)
            .map_err(|e| format!("Could not write genome manifest '{}': {e}", path.display()))
    }

    fn validate_manifest_files(manifest: &GenomeInstallManifest) -> Result<(), String> {
        for path in [
            &manifest.sequence_path,
            &manifest.annotation_path,
            &manifest.fasta_index_path,
        ] {
            if !Path::new(path).exists() {
                return Err(format!(
                    "Genome installation is incomplete; missing file '{}'",
                    path
                ));
            }
        }
        Ok(())
    }
}

fn now_unix_ms() -> u128 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

fn canonical_or_display(path: &Path) -> String {
    fs::canonicalize(path)
        .unwrap_or_else(|_| path.to_path_buf())
        .to_string_lossy()
        .into_owned()
}

fn sanitize_for_path(s: &str) -> String {
    let mut out = String::new();
    for c in s.chars() {
        if c.is_ascii_alphanumeric() {
            out.push(c.to_ascii_lowercase());
        } else if matches!(c, ' ' | '-' | '_' | '.') {
            if !out.ends_with('_') {
                out.push('_');
            }
        }
    }
    let trimmed = out.trim_matches('_');
    if trimmed.is_empty() {
        "genome".to_string()
    } else {
        trimmed.to_string()
    }
}

fn resolve_ncbi_assembly_source(
    kind: &str,
    entry: &GenomeCatalogEntry,
) -> Result<Option<String>, String> {
    let accession = entry
        .ncbi_assembly_accession
        .as_ref()
        .map(|v| v.trim())
        .filter(|v| !v.is_empty());
    let assembly_name = entry
        .ncbi_assembly_name
        .as_ref()
        .map(|v| v.trim())
        .filter(|v| !v.is_empty());
    match (accession, assembly_name) {
        (None, None) => Ok(None),
        (Some(_), None) | (None, Some(_)) => Err(
            "NCBI assembly source requires both 'ncbi_assembly_accession' and 'ncbi_assembly_name'"
                .to_string(),
        ),
        (Some(accession), Some(assembly_name)) => {
            let normalized_accession = accession.trim().to_ascii_uppercase();
            let normalized_assembly_name = normalize_ncbi_assembly_name(assembly_name);
            let base =
                build_ncbi_assembly_ftp_base(&normalized_accession, &normalized_assembly_name)?;
            let suffix = match kind {
                "sequence" => "_genomic.fna.gz",
                "annotation" => "_genomic.gff.gz",
                _ => return Err(format!("Unsupported NCBI assembly source kind '{kind}'")),
            };
            Ok(Some(format!(
                "{base}/{}_{}{}",
                normalized_accession, normalized_assembly_name, suffix
            )))
        }
    }
}

fn resolve_genbank_accession_source(
    kind: &str,
    entry: &GenomeCatalogEntry,
) -> Result<Option<String>, String> {
    let accession = entry
        .genbank_accession
        .as_ref()
        .map(|v| v.trim())
        .filter(|v| !v.is_empty());
    let Some(raw_accession) = accession else {
        return Ok(None);
    };
    let accession = validate_genbank_accession(raw_accession)?;
    let rettype = match kind {
        "sequence" => "fasta",
        // Keep complete features/parts for vector-oriented annotation searches.
        "annotation" => "gbwithparts",
        _ => {
            return Err(format!(
                "Unsupported GenBank accession source kind '{kind}'"
            ))
        }
    };
    Ok(Some(format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype={rettype}&retmode=text"
    )))
}

fn validate_genbank_accession(raw: &str) -> Result<String, String> {
    let value = raw.trim();
    if value.is_empty() {
        return Err("GenBank accession is empty".to_string());
    }
    let upper = value.to_ascii_uppercase();
    if matches!(
        upper.as_str(),
        "LOCAL_UNPUBLISHED" | "NOT_UPLOADED_TO_GENBANK" | "NOT_UPLOADED"
    ) {
        return Err(format!(
            "GenBank accession '{value}' is marked as unpublished/local-only and cannot be fetched remotely"
        ));
    }
    let first = value.chars().next().unwrap_or_default();
    if !first.is_ascii_alphabetic() {
        return Err(format!(
            "Invalid GenBank accession '{value}' (must start with an ASCII letter)"
        ));
    }
    if !value
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-'))
    {
        return Err(format!(
            "Invalid GenBank accession '{value}' (allowed characters: A-Z, 0-9, '_', '.', '-')"
        ));
    }
    Ok(value.to_string())
}

fn build_ncbi_assembly_ftp_base(accession: &str, assembly_name: &str) -> Result<String, String> {
    let normalized_accession = accession.trim().to_ascii_uppercase();
    let Some((prefix, rest)) = normalized_accession.split_once('_') else {
        return Err(format!(
            "Invalid NCBI assembly accession '{accession}' (expected GCA_/GCF_ prefix)"
        ));
    };
    if prefix != "GCA" && prefix != "GCF" {
        return Err(format!(
            "Invalid NCBI assembly accession '{accession}' (expected GCA_/GCF_ prefix)"
        ));
    }
    let Some((digits, version)) = rest.split_once('.') else {
        return Err(format!(
            "Invalid NCBI assembly accession '{accession}' (expected version suffix like '.2')"
        ));
    };
    if digits.is_empty()
        || !digits.chars().all(|c| c.is_ascii_digit())
        || version.is_empty()
        || !version.chars().all(|c| c.is_ascii_digit())
    {
        return Err(format!(
            "Invalid NCBI assembly accession '{accession}' (non-numeric accession/version segment)"
        ));
    }
    let mut padded_digits = digits.to_string();
    let rem = padded_digits.len() % 3;
    if rem != 0 {
        padded_digits = format!("{}{}", "0".repeat(3 - rem), padded_digits);
    }
    let segments: Vec<String> = padded_digits
        .as_bytes()
        .chunks(3)
        .map(|chunk| String::from_utf8_lossy(chunk).to_string())
        .collect();
    let assembly_name = normalize_ncbi_assembly_name(assembly_name);
    if assembly_name.is_empty() || assembly_name.contains('/') || assembly_name.contains('\\') {
        return Err(format!(
            "Invalid NCBI assembly name '{assembly_name}' for accession '{accession}'"
        ));
    }
    Ok(format!(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}_{}",
        prefix,
        segments.join("/"),
        normalized_accession,
        assembly_name
    ))
}

fn normalize_ncbi_assembly_name(assembly_name: &str) -> String {
    assembly_name
        .trim()
        .chars()
        .map(|c| if c.is_whitespace() { '_' } else { c })
        .collect::<String>()
}

fn infer_annotation_extension(source: &str) -> &'static str {
    let lower = source.to_ascii_lowercase();
    if lower.contains(".gbff")
        || lower.contains(".gbk")
        || lower.contains(".genbank")
        || lower.contains(".gb")
        || lower.contains("rettype=gbwithparts")
        || lower.contains("rettype=gb")
        || lower.contains("rettype=genbank")
    {
        "gbff"
    } else if lower.contains(".gff3") {
        "gff3"
    } else if lower.contains(".gff") {
        "gff"
    } else {
        "gtf"
    }
}

fn is_http_source(source: &str) -> bool {
    source.starts_with("http://") || source.starts_with("https://")
}

fn is_gzip_source(source: &str) -> bool {
    source.to_ascii_lowercase().ends_with(".gz")
}

struct SourceReader {
    reader: Box<dyn Read>,
    total_bytes: Option<u64>,
}

fn open_source_reader(source: &str) -> Result<SourceReader, String> {
    if is_http_source(source) {
        let response = get(source)
            .map_err(|e| format!("Could not fetch '{source}': {e}"))?
            .error_for_status()
            .map_err(|e| format!("Could not fetch '{source}': {e}"))?;
        let total_bytes = response.content_length();
        return Ok(SourceReader {
            reader: Box::new(response),
            total_bytes,
        });
    }
    let path = if let Some(stripped) = source.strip_prefix("file://") {
        PathBuf::from(stripped)
    } else {
        PathBuf::from(source)
    };
    let file = File::open(&path)
        .map_err(|e| format!("Could not open source file '{}': {e}", path.display()))?;
    let total_bytes = file.metadata().ok().map(|m| m.len());
    Ok(SourceReader {
        reader: Box::new(file),
        total_bytes,
    })
}

struct ProgressReader<R, F> {
    inner: R,
    callback: F,
    bytes_done: u64,
}

impl<R, F> ProgressReader<R, F> {
    fn new(inner: R, callback: F) -> Self {
        Self {
            inner,
            callback,
            bytes_done: 0,
        }
    }
}

impl<R: Read, F: FnMut(u64)> Read for ProgressReader<R, F> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.bytes_done += n as u64;
        (self.callback)(self.bytes_done);
        Ok(n)
    }
}

fn materialize_source_with_progress<F>(
    source: &str,
    destination: &Path,
    mut on_progress: F,
) -> Result<(), String>
where
    F: FnMut(u64, Option<u64>),
{
    if let Some(parent) = destination.parent() {
        fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Could not create destination directory '{}': {e}",
                parent.display()
            )
        })?;
    }
    let mut tmp_os: OsString = destination.as_os_str().to_os_string();
    tmp_os.push(".part");
    let tmp_path = PathBuf::from(tmp_os);

    let SourceReader {
        reader,
        total_bytes,
    } = open_source_reader(source)?;
    on_progress(0, total_bytes);
    let mut writer = BufWriter::new(
        File::create(&tmp_path)
            .map_err(|e| format!("Could not create '{}': {e}", tmp_path.display()))?,
    );

    let copy_result = if is_gzip_source(source) {
        let progress_reader = ProgressReader::new(reader, |done| on_progress(done, total_bytes));
        let mut decoder = GzDecoder::new(progress_reader);
        std::io::copy(&mut decoder, &mut writer)
            .map_err(|e| format!("Could not decompress '{source}': {e}"))
    } else {
        let mut progress_reader =
            ProgressReader::new(reader, |done| on_progress(done, total_bytes));
        std::io::copy(&mut progress_reader, &mut writer)
            .map_err(|e| format!("Could not copy '{source}': {e}"))
    };

    if let Err(e) = copy_result {
        let _ = fs::remove_file(&tmp_path);
        return Err(e);
    }

    writer
        .flush()
        .map_err(|e| format!("Could not flush '{}': {e}", tmp_path.display()))?;
    let final_done = fs::metadata(&tmp_path).ok().map(|m| m.len()).unwrap_or(0);
    on_progress(total_bytes.unwrap_or(final_done), total_bytes);
    fs::rename(&tmp_path, destination).map_err(|e| {
        format!(
            "Could not finalize destination '{}': {e}",
            destination.display()
        )
    })?;
    Ok(())
}

fn build_fasta_index_with_progress<F>(
    fasta_path: &Path,
    index_path: &Path,
    mut on_progress: F,
) -> Result<(), String>
where
    F: FnMut(u64, Option<u64>),
{
    let total_bytes = fs::metadata(fasta_path).ok().map(|m| m.len());
    let file = File::open(fasta_path)
        .map_err(|e| format!("Could not open FASTA '{}': {e}", fasta_path.display()))?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let mut byte_offset: u64 = 0;
    let mut entries: Vec<(String, FastaIndexEntry, bool)> = Vec::new();
    // tuple fields: (name, entry, saw_short_line)
    let mut active: Option<(String, FastaIndexEntry, bool)> = None;
    on_progress(0, total_bytes);

    loop {
        line.clear();
        let bytes_read = reader
            .read_line(&mut line)
            .map_err(|e| format!("Could not read FASTA '{}': {e}", fasta_path.display()))?;
        if bytes_read == 0 {
            break;
        }
        let stripped = line.trim_end_matches(&['\n', '\r'][..]);

        if stripped.starts_with('>') {
            if let Some((name, entry, _)) = active.take() {
                if entry.length == 0 {
                    return Err(format!(
                        "FASTA '{}' has empty sequence record '{}'",
                        fasta_path.display(),
                        name
                    ));
                }
                entries.push((name, entry, false));
            }
            let name = stripped
                .trim_start_matches('>')
                .split_whitespace()
                .next()
                .ok_or_else(|| {
                    format!("FASTA '{}' has malformed header line", fasta_path.display())
                })?
                .to_string();
            active = Some((
                name,
                FastaIndexEntry {
                    length: 0,
                    offset: 0,
                    line_bases: 0,
                    line_bytes: 0,
                },
                false,
            ));
        } else if !stripped.is_empty() {
            let seq_bases = stripped.len() as u64;
            if let Some((_, ref mut entry, ref mut saw_short_line)) = active {
                if entry.line_bases == 0 {
                    entry.offset = byte_offset;
                    entry.line_bases = seq_bases;
                    entry.line_bytes = bytes_read as u64;
                } else {
                    if seq_bases > entry.line_bases {
                        return Err(format!(
                            "FASTA '{}' has inconsistent line length; line longer than first line",
                            fasta_path.display()
                        ));
                    }
                    if *saw_short_line && seq_bases == entry.line_bases {
                        return Err(format!(
                            "FASTA '{}' has inconsistent line length after a short line",
                            fasta_path.display()
                        ));
                    }
                    if seq_bases < entry.line_bases {
                        *saw_short_line = true;
                    }
                }
                entry.length += seq_bases;
            } else {
                return Err(format!(
                    "FASTA '{}' contains sequence data before first header",
                    fasta_path.display()
                ));
            }
        }
        byte_offset += bytes_read as u64;
        on_progress(byte_offset, total_bytes);
    }

    if let Some((name, entry, _)) = active {
        if entry.length == 0 {
            return Err(format!(
                "FASTA '{}' has empty sequence record '{}'",
                fasta_path.display(),
                name
            ));
        }
        entries.push((name, entry, false));
    }
    if entries.is_empty() {
        return Err(format!(
            "FASTA '{}' does not contain any sequence records",
            fasta_path.display()
        ));
    }

    if let Some(parent) = index_path.parent() {
        fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Could not create FASTA index directory '{}': {e}",
                parent.display()
            )
        })?;
    }
    let mut writer = BufWriter::new(File::create(index_path).map_err(|e| {
        format!(
            "Could not create FASTA index '{}': {e}",
            index_path.display()
        )
    })?);
    for (name, entry, _) in entries {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            name, entry.length, entry.offset, entry.line_bases, entry.line_bytes
        )
        .map_err(|e| {
            format!(
                "Could not write FASTA index '{}': {e}",
                index_path.display()
            )
        })?;
    }
    writer.flush().map_err(|e| {
        format!(
            "Could not flush FASTA index '{}': {e}",
            index_path.display()
        )
    })?;
    on_progress(total_bytes.unwrap_or(byte_offset), total_bytes);
    Ok(())
}

fn load_fasta_index(path: &Path) -> Result<HashMap<String, FastaIndexEntry>, String> {
    let file = File::open(path)
        .map_err(|e| format!("Could not open FASTA index '{}': {e}", path.display()))?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    for (i, line) in reader.lines().enumerate() {
        let line =
            line.map_err(|e| format!("Could not read FASTA index '{}': {e}", path.display()))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 5 {
            return Err(format!(
                "Invalid FASTA index line {} in '{}': expected 5 tab-separated fields",
                i + 1,
                path.display()
            ));
        }
        let name = cols[0].to_string();
        let length = cols[1].parse::<u64>().map_err(|e| {
            format!(
                "Invalid FASTA index length '{}' at line {} in '{}': {e}",
                cols[1],
                i + 1,
                path.display()
            )
        })?;
        let offset = cols[2].parse::<u64>().map_err(|e| {
            format!(
                "Invalid FASTA index offset '{}' at line {} in '{}': {e}",
                cols[2],
                i + 1,
                path.display()
            )
        })?;
        let line_bases = cols[3].parse::<u64>().map_err(|e| {
            format!(
                "Invalid FASTA index line_bases '{}' at line {} in '{}': {e}",
                cols[3],
                i + 1,
                path.display()
            )
        })?;
        let line_bytes = cols[4].parse::<u64>().map_err(|e| {
            format!(
                "Invalid FASTA index line_bytes '{}' at line {} in '{}': {e}",
                cols[4],
                i + 1,
                path.display()
            )
        })?;
        if line_bases == 0 || line_bytes == 0 {
            return Err(format!(
                "Invalid FASTA index line {} in '{}': line_bases/line_bytes must be > 0",
                i + 1,
                path.display()
            ));
        }
        map.insert(
            name,
            FastaIndexEntry {
                length,
                offset,
                line_bases,
                line_bytes,
            },
        );
    }
    if map.is_empty() {
        return Err(format!("FASTA index '{}' is empty", path.display()));
    }
    Ok(map)
}

fn read_fasta_slice(
    fasta_path: &Path,
    entry: &FastaIndexEntry,
    start_1based: u64,
    end_1based: u64,
) -> Result<String, String> {
    let start0 = start_1based - 1;
    let target_len = (end_1based - start_1based + 1) as usize;
    let row = start0 / entry.line_bases;
    let col = start0 % entry.line_bases;
    let seek_pos = entry
        .offset
        .checked_add(row.saturating_mul(entry.line_bytes))
        .and_then(|v| v.checked_add(col))
        .ok_or_else(|| "Computed FASTA seek position overflowed".to_string())?;

    let mut file = File::open(fasta_path)
        .map_err(|e| format!("Could not open FASTA '{}': {e}", fasta_path.display()))?;
    file.seek(SeekFrom::Start(seek_pos))
        .map_err(|e| format!("Could not seek FASTA '{}': {e}", fasta_path.display()))?;
    let mut reader = BufReader::new(file);

    let mut out: Vec<u8> = Vec::with_capacity(target_len);
    let mut chunk = [0u8; 8192];
    while out.len() < target_len {
        let n = reader
            .read(&mut chunk)
            .map_err(|e| format!("Could not read FASTA '{}': {e}", fasta_path.display()))?;
        if n == 0 {
            break;
        }
        for b in &chunk[..n] {
            if *b == b'\n' || *b == b'\r' {
                continue;
            }
            out.push(*b);
            if out.len() == target_len {
                break;
            }
        }
    }
    if out.len() != target_len {
        return Err(format!(
            "Could not read requested interval from FASTA '{}'; expected {} bases, got {}",
            fasta_path.display(),
            target_len,
            out.len()
        ));
    }
    String::from_utf8(out).map_err(|e| format!("Extracted sequence is not valid UTF-8: {e}"))
}

fn build_gene_index_file<F>(
    annotation_path: &Path,
    gene_index_path: &Path,
    on_progress: F,
) -> Result<(), String>
where
    F: FnMut(u64, Option<u64>),
{
    let genes = parse_annotation_gene_records_with_progress(annotation_path, on_progress)?;
    if let Some(parent) = gene_index_path.parent() {
        fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Could not create gene-index directory '{}': {e}",
                parent.display()
            )
        })?;
    }
    let file = File::create(gene_index_path).map_err(|e| {
        format!(
            "Could not create gene-index file '{}': {e}",
            gene_index_path.display()
        )
    })?;
    serde_json::to_writer(BufWriter::new(file), &genes).map_err(|e| {
        format!(
            "Could not write gene index '{}': {e}",
            gene_index_path.display()
        )
    })?;
    Ok(())
}

fn load_gene_index_file(gene_index_path: &Path) -> Result<Vec<GenomeGeneRecord>, String> {
    let file = File::open(gene_index_path).map_err(|e| {
        format!(
            "Could not open gene index '{}': {e}",
            gene_index_path.display()
        )
    })?;
    serde_json::from_reader(BufReader::new(file)).map_err(|e| {
        format!(
            "Could not parse gene index '{}': {e}",
            gene_index_path.display()
        )
    })
}

fn parse_annotation_gene_records_with_progress<F>(
    path: &Path,
    on_progress: F,
) -> Result<Vec<GenomeGeneRecord>, String>
where
    F: FnMut(u64, Option<u64>),
{
    if is_genbank_annotation_path(path) {
        parse_genbank_gene_records_with_progress(path, on_progress)
    } else {
        parse_tabular_annotation_gene_records_with_progress(path, on_progress)
    }
}

fn is_genbank_annotation_path(path: &Path) -> bool {
    let lower = path
        .file_name()
        .and_then(|v| v.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    lower.ends_with(".gb")
        || lower.ends_with(".gbk")
        || lower.ends_with(".gbff")
        || lower.ends_with(".genbank")
}

fn parse_tabular_annotation_gene_records_with_progress<F>(
    path: &Path,
    mut on_progress: F,
) -> Result<Vec<GenomeGeneRecord>, String>
where
    F: FnMut(u64, Option<u64>),
{
    let total_bytes = fs::metadata(path).ok().map(|m| m.len());
    let file = File::open(path)
        .map_err(|e| format!("Could not open annotation file '{}': {e}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut genes: Vec<GenomeGeneRecord> = vec![];
    let mut bytes_read_total: u64 = 0;
    let mut raw_line: Vec<u8> = vec![];
    on_progress(0, total_bytes);
    loop {
        raw_line.clear();
        let line_bytes = reader
            .read_until(b'\n', &mut raw_line)
            .map_err(|e| format!("Could not read annotation file '{}': {e}", path.display()))?;
        if line_bytes == 0 {
            break;
        }
        bytes_read_total = bytes_read_total.saturating_add(line_bytes as u64);
        on_progress(bytes_read_total, total_bytes);
        if raw_line.ends_with(b"\n") {
            raw_line.pop();
        }
        if raw_line.ends_with(b"\r") {
            raw_line.pop();
        }
        let line = String::from_utf8_lossy(&raw_line);
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let mut cols = trimmed.splitn(9, '\t');
        let chromosome = match cols.next().map(str::trim) {
            Some(v) if !v.is_empty() => v,
            _ => continue,
        };
        if cols.next().is_none() {
            continue;
        }
        let feature_kind = match cols.next() {
            Some(v) => v,
            None => continue,
        };
        let start_1based = match cols.next().and_then(parse_annotation_coordinate) {
            Some(v) => v,
            None => continue,
        };
        let end_1based = match cols.next().and_then(parse_annotation_coordinate) {
            Some(v) => v,
            None => continue,
        };
        if cols.next().is_none() {
            continue;
        }
        let strand_raw = match cols.next() {
            Some(v) => v,
            None => continue,
        };
        if cols.next().is_none() {
            continue;
        }
        let attrs_raw = match cols.next() {
            Some(v) => v,
            None => continue,
        };

        if !is_gene_feature_kind(feature_kind) {
            continue;
        }
        if start_1based == 0 || end_1based < start_1based {
            continue;
        }
        let attrs = parse_annotation_attributes(attrs_raw);
        let gene_name = pick_gene_name(&attrs);
        let mut gene_id = pick_gene_id(&attrs);
        let biotype = pick_gene_biotype(&attrs);
        if gene_name.is_none() && gene_id.is_none() {
            gene_id = Some(format!("{}_{}_{}", chromosome, start_1based, end_1based));
        }
        let strand = strand_raw.chars().next().and_then(|c| match c {
            '+' | '-' => Some(c),
            _ => None,
        });
        genes.push(GenomeGeneRecord {
            chromosome: chromosome.to_string(),
            start_1based,
            end_1based,
            strand,
            gene_id,
            gene_name,
            biotype,
        });
    }
    genes.sort_by(|a, b| {
        let a_name = a
            .gene_name
            .as_deref()
            .or(a.gene_id.as_deref())
            .unwrap_or(&a.chromosome);
        let b_name = b
            .gene_name
            .as_deref()
            .or(b.gene_id.as_deref())
            .unwrap_or(&b.chromosome);
        a_name
            .cmp(b_name)
            .then(a.chromosome.cmp(&b.chromosome))
            .then(a.start_1based.cmp(&b.start_1based))
            .then(a.end_1based.cmp(&b.end_1based))
    });
    on_progress(total_bytes.unwrap_or(bytes_read_total), total_bytes);
    Ok(genes)
}

fn parse_genbank_gene_records_with_progress<F>(
    path: &Path,
    mut on_progress: F,
) -> Result<Vec<GenomeGeneRecord>, String>
where
    F: FnMut(u64, Option<u64>),
{
    let total_bytes = fs::metadata(path).ok().map(|m| m.len());
    on_progress(0, total_bytes);
    let parsed = gb_io::reader::parse_file(path.to_string_lossy().as_ref()).map_err(|e| {
        format!(
            "Could not parse GenBank annotation file '{}': {e}",
            path.display()
        )
    })?;

    let mut records: Vec<GenomeGeneRecord> = vec![];
    for (record_idx, seq) in parsed.iter().enumerate() {
        let chromosome = genbank_record_chromosome(seq, record_idx);
        for feature in &seq.features {
            let kind = feature.kind.to_string();
            if !is_genbank_feature_kind_indexed(&kind) {
                continue;
            }
            let Ok((raw_from, raw_to)) = feature.location.find_bounds() else {
                continue;
            };
            if raw_from < 0 || raw_to < 0 {
                continue;
            }
            let mut start_0based = raw_from as usize;
            let mut end_0based = raw_to as usize;
            if end_0based < start_0based {
                std::mem::swap(&mut start_0based, &mut end_0based);
            }
            let attrs = genbank_feature_attributes(feature);
            let biotype = normalize_genbank_feature_biotype(&kind, &attrs);
            let gene_name = pick_gene_name(&attrs).or_else(|| {
                pick_annotation_attribute(
                    &attrs,
                    &[
                        "label",
                        "product",
                        "note",
                        "standard_name",
                        "bound_moiety",
                        "regulatory_class",
                    ],
                )
            });
            let mut gene_id = pick_gene_id(&attrs);
            if gene_id.is_none() {
                gene_id = Some(format!(
                    "{}:{}:{}-{}",
                    kind.to_ascii_lowercase(),
                    chromosome,
                    start_0based + 1,
                    end_0based + 1
                ));
            }
            let strand = Some(if feature_is_reverse(feature) {
                '-'
            } else {
                '+'
            });
            records.push(GenomeGeneRecord {
                chromosome: chromosome.clone(),
                start_1based: start_0based + 1,
                end_1based: end_0based + 1,
                strand,
                gene_id,
                gene_name,
                biotype,
            });
        }
    }
    records.sort_by(|a, b| {
        let a_name = a
            .gene_name
            .as_deref()
            .or(a.gene_id.as_deref())
            .unwrap_or(&a.chromosome);
        let b_name = b
            .gene_name
            .as_deref()
            .or(b.gene_id.as_deref())
            .unwrap_or(&b.chromosome);
        a_name
            .cmp(b_name)
            .then(a.chromosome.cmp(&b.chromosome))
            .then(a.start_1based.cmp(&b.start_1based))
            .then(a.end_1based.cmp(&b.end_1based))
    });
    let total_done = total_bytes.unwrap_or(0);
    on_progress(total_done, total_bytes);
    Ok(records)
}

fn genbank_record_chromosome(seq: &gb_io::seq::Seq, record_idx: usize) -> String {
    let from_version = seq
        .version
        .as_ref()
        .map(|v| v.trim())
        .filter(|v| !v.is_empty())
        .map(str::to_string);
    let from_accession = seq
        .accession
        .as_ref()
        .map(|v| v.trim())
        .filter(|v| !v.is_empty())
        .map(str::to_string);
    let from_name = seq
        .name
        .as_ref()
        .map(|v| v.trim())
        .filter(|v| !v.is_empty())
        .map(str::to_string);
    from_version
        .or(from_accession)
        .or(from_name)
        .unwrap_or_else(|| format!("record_{}", record_idx + 1))
}

fn is_genbank_feature_kind_indexed(raw: &str) -> bool {
    let lower = raw.trim().to_ascii_lowercase();
    if lower.is_empty() {
        return false;
    }
    if lower == "source" {
        return false;
    }
    true
}

fn genbank_feature_attributes(feature: &gb_io::seq::Feature) -> HashMap<String, String> {
    let mut attrs = HashMap::new();
    for (key, value) in &feature.qualifiers {
        let key = key.trim().to_ascii_lowercase();
        if key.is_empty() {
            continue;
        }
        let Some(raw) = value.as_ref() else {
            continue;
        };
        if let Some(normalized) = normalize_annotation_attribute_value(raw) {
            attrs.entry(key).or_insert(normalized);
        }
    }
    attrs
}

fn normalize_genbank_feature_biotype(
    feature_kind: &str,
    attrs: &HashMap<String, String>,
) -> Option<String> {
    let kind = feature_kind.trim().to_ascii_lowercase();
    if kind.is_empty() {
        return None;
    }
    if kind == "regulatory" {
        if let Some(class) = attrs.get("regulatory_class").map(|v| v.trim()) {
            if !class.is_empty() {
                return Some(class.to_ascii_lowercase());
            }
        }
    }
    if kind == "rep_origin" {
        return Some("origin_of_replication".to_string());
    }
    if kind == "misc_feature" {
        for key in ["label", "note"] {
            if let Some(value) = attrs.get(key) {
                let lower = value.to_ascii_lowercase();
                if lower.contains("ltr") {
                    return Some("ltr".to_string());
                }
                if lower.contains("itr") {
                    return Some("itr".to_string());
                }
                if lower.contains("origin") || lower.contains(" ori") || lower.starts_with("ori") {
                    return Some("origin_of_replication".to_string());
                }
                if lower.contains("promoter") {
                    return Some("promoter".to_string());
                }
                if lower.contains("marker") || lower.contains("selection") {
                    return Some("selectable_marker".to_string());
                }
            }
        }
    }
    Some(kind)
}

fn is_gene_feature_kind(raw: &str) -> bool {
    let lower = raw.trim().to_ascii_lowercase();
    lower == "gene" || lower == "pseudogene" || lower.ends_with("_gene")
}

fn parse_annotation_coordinate(raw: &str) -> Option<usize> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    let stripped = trimmed.trim_start_matches(['<', '>', '=']);
    if stripped.is_empty() {
        return None;
    }
    let no_commas = stripped.replace(',', "");
    if let Ok(v) = no_commas.parse::<usize>() {
        return (v > 0).then_some(v);
    }
    let prefix_digits: String = no_commas
        .chars()
        .take_while(|ch| ch.is_ascii_digit())
        .collect();
    if prefix_digits.is_empty() {
        None
    } else {
        prefix_digits.parse::<usize>().ok().filter(|v| *v > 0)
    }
}

fn pick_gene_name(attrs: &HashMap<String, String>) -> Option<String> {
    pick_annotation_attribute(
        attrs,
        &[
            "gene_name",
            "name",
            "gene_symbol",
            "symbol",
            "gene",
            "locus_tag",
            "standard_name",
        ],
    )
}

fn pick_gene_id(attrs: &HashMap<String, String>) -> Option<String> {
    let direct = pick_annotation_attribute(
        attrs,
        &[
            "gene_id",
            "id",
            "gene",
            "locus_tag",
            "geneid",
            "transcript_id",
        ],
    );
    direct
        .or_else(|| pick_gene_id_from_dbxref(attrs))
        .map(normalize_gene_id)
}

fn pick_gene_biotype(attrs: &HashMap<String, String>) -> Option<String> {
    pick_annotation_attribute(
        attrs,
        &[
            "gene_biotype",
            "gene_type",
            "biotype",
            "transcript_biotype",
            "transcript_type",
        ],
    )
}

fn pick_gene_id_from_dbxref(attrs: &HashMap<String, String>) -> Option<String> {
    let raw = attrs
        .get("dbxref")
        .or_else(|| attrs.get("db_xref"))
        .or_else(|| attrs.get("xref"))?;
    for token in raw.split(',').map(str::trim).filter(|v| !v.is_empty()) {
        if let Some((key, value)) = token.split_once(':') {
            let k = key.trim().to_ascii_lowercase();
            let v = value.trim();
            if v.is_empty() {
                continue;
            }
            if matches!(
                k.as_str(),
                "geneid" | "gene" | "ensembl_gene_id" | "ensemblgene" | "hgnc" | "locus_tag"
            ) {
                return Some(v.to_string());
            }
        } else {
            return Some(token.to_string());
        }
    }
    None
}

fn normalize_gene_id(raw: String) -> String {
    let trimmed = raw.trim();
    if let Some((prefix, rest)) = trimmed.split_once(':') {
        let p = prefix.trim().to_ascii_lowercase();
        if matches!(p.as_str(), "gene" | "geneid" | "id") {
            let rest = rest.trim();
            if !rest.is_empty() {
                return rest.to_string();
            }
        }
    }
    trimmed.to_string()
}

fn pick_annotation_attribute(attrs: &HashMap<String, String>, keys: &[&str]) -> Option<String> {
    keys.iter()
        .filter_map(|key| attrs.get(*key))
        .find(|value| !value.trim().is_empty())
        .cloned()
}

fn split_annotation_attributes(raw: &str) -> Vec<String> {
    let mut fields: Vec<String> = vec![];
    let mut current = String::new();
    let mut in_quotes = false;
    let mut escape = false;
    for ch in raw.chars() {
        if escape {
            current.push(ch);
            escape = false;
            continue;
        }
        match ch {
            '\\' if in_quotes => {
                current.push(ch);
                escape = true;
            }
            '"' => {
                in_quotes = !in_quotes;
                current.push(ch);
            }
            ';' if !in_quotes => {
                let field = current.trim();
                if !field.is_empty() {
                    fields.push(field.to_string());
                }
                current.clear();
            }
            _ => current.push(ch),
        }
    }
    let tail = current.trim();
    if !tail.is_empty() {
        fields.push(tail.to_string());
    }
    fields
}

fn normalize_annotation_attribute_value(raw: &str) -> Option<String> {
    let mut cleaned = raw
        .trim()
        .trim_matches('"')
        .trim()
        .trim_end_matches(',')
        .trim()
        .to_string();
    cleaned = cleaned
        .replace("\\\"", "\"")
        .replace("\\;", ";")
        .replace("\\,", ",");
    cleaned = cleaned.trim().trim_matches('"').trim().to_string();
    if cleaned.is_empty() {
        return None;
    }
    if cleaned.is_empty() {
        None
    } else {
        Some(cleaned)
    }
}

fn parse_annotation_attributes(raw: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for part in split_annotation_attributes(raw) {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if let Some((key, value)) = part.split_once('=') {
            let key = key.trim().to_ascii_lowercase();
            if !key.is_empty() {
                if let Some(value) = normalize_annotation_attribute_value(value) {
                    map.entry(key).or_insert(value);
                }
            }
            continue;
        }
        let mut pieces = part.splitn(2, char::is_whitespace);
        let key = pieces
            .next()
            .unwrap_or_default()
            .trim()
            .to_ascii_lowercase();
        let value = pieces.next().unwrap_or_default();
        if !key.is_empty() {
            if let Some(value) = normalize_annotation_attribute_value(value) {
                map.entry(key).or_insert(value);
            }
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{write::GzEncoder, Compression};
    use tempfile::tempdir;

    fn write_gzip(path: &Path, text: &str) {
        let file = File::create(path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(text.as_bytes()).unwrap();
        encoder.finish().unwrap();
    }

    fn file_url(path: &Path) -> String {
        format!("file://{}", path.display())
    }

    #[test]
    fn test_prepare_genome_once_and_extract_region() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta_gz = root.join("toy.fa.gz");
        let ann_gz = root.join("toy.gtf.gz");
        write_gzip(&fasta_gz, ">chr1\nACGT\nACGT\nACGT\n>II\nTTTT\nGGGG\n");
        write_gzip(
            &ann_gz,
            "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; gene_biotype \"protein_coding\";\n",
        );

        let cache_dir = root.join("cache");
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
    "sequence_remote": "{}",
    "annotations_remote": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            file_url(&fasta_gz),
            file_url(&ann_gz),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();

        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let first = catalog.prepare_genome_once("ToyGenome").unwrap();
        assert!(!first.reused_existing);
        assert!(Path::new(&first.sequence_path).exists());
        assert!(Path::new(&first.annotation_path).exists());
        assert!(Path::new(&first.fasta_index_path).exists());
        let manifest =
            GenomeCatalog::load_manifest(&cache_dir.join("toygenome").join("manifest.json"))
                .unwrap();
        let gene_index_path = manifest
            .gene_index_path
            .expect("gene index path should be set");
        assert!(Path::new(&gene_index_path).exists());

        let second = catalog.prepare_genome_once("ToyGenome").unwrap();
        assert!(second.reused_existing);

        let seq = catalog
            .get_sequence_region("ToyGenome", "1", 3, 10)
            .unwrap();
        assert_eq!(seq, "GTACGTAC");

        let genes = catalog.list_gene_regions("ToyGenome", None).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].chromosome, "chr1");
        assert_eq!(genes[0].start_1based, 1);
        assert_eq!(genes[0].end_1based, 12);
        assert_eq!(genes[0].strand, Some('+'));
        assert_eq!(genes[0].gene_id.as_deref(), Some("GENE1"));
        assert_eq!(genes[0].gene_name.as_deref(), Some("MYGENE"));
        assert_eq!(genes[0].biotype.as_deref(), Some("protein_coding"));
    }

    #[test]
    fn test_prepare_fails_when_annotation_source_missing() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta_gz = root.join("toy.fa.gz");
        write_gzip(&fasta_gz, ">chr1\nACGT\n");

        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
    "sequence_remote": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            file_url(&fasta_gz),
            root.join("cache").display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();

        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let err = catalog.prepare_genome_once("ToyGenome").unwrap_err();
        assert!(err.contains("annotation"));
    }

    #[test]
    fn test_parse_annotation_attributes_handles_quoted_semicolons() {
        let attrs = parse_annotation_attributes(
            r#"gene_id "GENE1"; gene_name "A;B"; Dbxref=GeneID:12345,HGNC:HGNC:5;"#,
        );
        assert_eq!(attrs.get("gene_id").map(String::as_str), Some("GENE1"));
        assert_eq!(attrs.get("gene_name").map(String::as_str), Some("A;B"));
        assert_eq!(
            attrs.get("dbxref").map(String::as_str),
            Some("GeneID:12345,HGNC:HGNC:5")
        );
    }

    #[test]
    fn test_parse_annotation_gene_records_tolerates_unexpected_lines() {
        let td = tempdir().unwrap();
        let path = td.path().join("mixed.gff");
        let mut bytes = vec![];
        bytes.extend_from_slice(b"##gff-version 3\n");
        bytes.extend_from_slice(
            b"chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene:GENE1;Name=Gene One;biotype=protein_coding\n",
        );
        bytes.extend_from_slice(
            b"chr1\tsrc\tprotein_coding_gene\t200\t300\t.\t-\t.\tgene_id \"GENE2\"; gene_name \"A;B\"; gene_biotype \"lncRNA\";\n",
        );
        bytes.extend_from_slice(b"chr1\tsrc\tgene\tbad\t400\t.\t+\t.\tID=gene:BAD\n");
        bytes.extend_from_slice(b"chr1\tsrc\tgene\t500\t400\t.\t+\t.\tID=gene:INV\n");
        bytes.extend_from_slice(
            b"chr1\tsrc\tpseudogene\t800\t900\t.\t+\t.\tDbxref=GeneID:12345,HGNC:HGNC:5\n",
        );
        bytes.extend_from_slice(b"chr1\tsrc\tgene\t<1001\t>1050\t.\t+\t.\tID=gene:GENE4\n");
        bytes.extend_from_slice(b"chr1\tsrc\tgene\t1200\t1300\t.\t+\t.\t.\n");
        bytes.extend_from_slice(&[0xff, 0xfe, 0x00, b'\n']);
        fs::write(&path, bytes).unwrap();

        let genes = parse_annotation_gene_records_with_progress(&path, |_done, _total| {}).unwrap();
        assert_eq!(genes.len(), 5);
        assert!(genes.iter().any(|g| g.gene_id.as_deref() == Some("GENE1")
            && g.gene_name.as_deref() == Some("Gene One")
            && g.biotype.as_deref() == Some("protein_coding")));
        assert!(genes.iter().any(|g| g.gene_id.as_deref() == Some("GENE2")
            && g.gene_name.as_deref() == Some("A;B")
            && g.biotype.as_deref() == Some("lncRNA")
            && g.strand == Some('-')));
        assert!(genes
            .iter()
            .any(|g| g.gene_id.as_deref() == Some("12345") && g.gene_name.is_none()));
        assert!(genes.iter().any(|g| g.gene_id.as_deref() == Some("GENE4")
            && g.start_1based == 1001
            && g.end_1based == 1050));
        assert!(genes
            .iter()
            .any(|g| g.gene_id.as_deref() == Some("chr1_1200_1300")));
    }

    #[test]
    fn test_build_ncbi_assembly_ftp_urls() {
        let entry = GenomeCatalogEntry {
            ncbi_assembly_accession: Some("gcf_000005845.2".to_string()),
            ncbi_assembly_name: Some("ASM584v2".to_string()),
            ..Default::default()
        };
        let seq = resolve_ncbi_assembly_source("sequence", &entry)
            .unwrap()
            .expect("sequence URL");
        let ann = resolve_ncbi_assembly_source("annotation", &entry)
            .unwrap()
            .expect("annotation URL");
        assert_eq!(
            seq,
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
        );
        assert_eq!(
            ann,
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz"
        );
    }

    #[test]
    fn test_ncbi_assembly_source_requires_accession_and_name() {
        let accession_only = GenomeCatalogEntry {
            ncbi_assembly_accession: Some("GCF_000005845.2".to_string()),
            ..Default::default()
        };
        let err = resolve_ncbi_assembly_source("sequence", &accession_only).unwrap_err();
        assert!(err.contains("requires both"));
    }

    #[test]
    fn test_build_genbank_efetch_urls() {
        let entry = GenomeCatalogEntry {
            genbank_accession: Some("L09137".to_string()),
            ..Default::default()
        };
        let seq = resolve_genbank_accession_source("sequence", &entry)
            .unwrap()
            .expect("sequence URL");
        let ann = resolve_genbank_accession_source("annotation", &entry)
            .unwrap()
            .expect("annotation URL");
        assert_eq!(
            seq,
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=L09137&rettype=fasta&retmode=text"
        );
        assert_eq!(
            ann,
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=L09137&rettype=gbwithparts&retmode=text"
        );
    }

    #[test]
    fn test_genbank_accession_rejects_unpublished_placeholders() {
        let entry = GenomeCatalogEntry {
            genbank_accession: Some("LOCAL_UNPUBLISHED".to_string()),
            ..Default::default()
        };
        let err = resolve_genbank_accession_source("sequence", &entry).unwrap_err();
        assert!(err.contains("unpublished"));
    }

    #[test]
    fn test_local_missing_falls_back_to_genbank_accession_source() {
        let entry = GenomeCatalogEntry {
            sequence_local: Some("missing/local.fa.gz".to_string()),
            genbank_accession: Some("L09137".to_string()),
            ..Default::default()
        };
        let mut entries = HashMap::new();
        entries.insert("Vec".to_string(), entry.clone());
        let catalog = GenomeCatalog {
            entries,
            catalog_base_dir: PathBuf::from("."),
        };
        let source = catalog
            .resolve_source(
                "Vec",
                "sequence",
                entry.sequence_local.as_ref(),
                entry.sequence_remote.as_ref(),
                &entry,
            )
            .unwrap();
        assert!(source.contains("db=nuccore&id=L09137"));
        assert!(source.contains("rettype=fasta"));
    }

    #[test]
    fn test_parse_genbank_annotation_records_pgex() {
        let genes = parse_annotation_gene_records_with_progress(
            Path::new("test_files/pGEX-3X.gb"),
            |_done, _total| {},
        )
        .unwrap();
        assert!(!genes.is_empty());
        assert!(genes
            .iter()
            .any(|g| g.biotype.as_deref() == Some("promoter")));
        assert!(genes.iter().any(|g| {
            g.biotype.as_deref() == Some("origin_of_replication")
                && g.start_1based <= 2978
                && g.end_1based >= 2978
        }));
        assert!(genes.iter().any(|g| {
            let is_bla = g
                .gene_name
                .as_ref()
                .map(|name| name.eq_ignore_ascii_case("bla"))
                .unwrap_or(false)
                || g.gene_id
                    .as_ref()
                    .map(|id| id.eq_ignore_ascii_case("bla"))
                    .unwrap_or(false);
            is_bla
                && (1289..=1291).contains(&g.start_1based)
                && (2219..=2221).contains(&g.end_1based)
        }));
    }

    #[test]
    fn test_normalize_genbank_biotypes_for_helper_vectors() {
        let mut attrs = HashMap::new();
        attrs.insert("label".to_string(), "5' LTR".to_string());
        assert_eq!(
            normalize_genbank_feature_biotype("misc_feature", &attrs).as_deref(),
            Some("ltr")
        );

        attrs.insert("label".to_string(), "AAV ITR left".to_string());
        assert_eq!(
            normalize_genbank_feature_biotype("misc_feature", &attrs).as_deref(),
            Some("itr")
        );

        attrs.insert("label".to_string(), "pUC ori".to_string());
        assert_eq!(
            normalize_genbank_feature_biotype("misc_feature", &attrs).as_deref(),
            Some("origin_of_replication")
        );

        attrs.insert("label".to_string(), "CMV promoter".to_string());
        assert_eq!(
            normalize_genbank_feature_biotype("misc_feature", &attrs).as_deref(),
            Some("promoter")
        );
    }
}
