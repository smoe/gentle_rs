use flate2::read::GzDecoder;
use reqwest::blocking::get;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

pub const DEFAULT_GENOME_CATALOG_PATH: &str = "assets/genomes.json";
pub const DEFAULT_GENOME_CACHE_DIR: &str = "data/genomes";

/// Catalog entry describing where to fetch one genome assembly and annotation.
#[derive(Default, Deserialize, Serialize, Debug, Clone)]
pub struct GenomeCatalogEntry {
    pub ncbi_taxonomy_id: Option<u32>,
    pub description: Option<String>,
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
        Ok(Self::validate_manifest_files(&manifest).is_ok())
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
            let manifest = Self::load_manifest(&manifest_path)?;
            Self::validate_manifest_files(&manifest)?;
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
        )?;
        let annotation_source = self.resolve_source(
            genome_id,
            "annotation",
            entry.annotations_local.as_ref(),
            entry.annotations_remote.as_ref(),
        )?;

        let sequence_path = install_dir.join("sequence.fa");
        let annotation_ext = infer_annotation_extension(&annotation_source);
        let annotation_path = install_dir.join(format!("annotation.{annotation_ext}"));
        let fasta_index_path = install_dir.join("sequence.fa.fai");

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

        let manifest = GenomeInstallManifest {
            genome_id: genome_id.to_string(),
            sequence_source,
            annotation_source,
            sequence_path: canonical_or_display(&sequence_path),
            annotation_path: canonical_or_display(&annotation_path),
            fasta_index_path: canonical_or_display(&fasta_index_path),
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
        parse_annotation_gene_records(Path::new(&manifest.annotation_path))
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
            .or_else(|| entry.cache_dir.as_ref().map(|raw| self.resolve_local_path(raw)))
            .unwrap_or_else(|| PathBuf::from(DEFAULT_GENOME_CACHE_DIR));
        base.join(sanitize_for_path(genome_id))
    }

    fn resolve_source(
        &self,
        genome_id: &str,
        kind: &str,
        local: Option<&String>,
        remote: Option<&String>,
    ) -> Result<String, String> {
        if let Some(local_raw) = local {
            let local_path = self.resolve_local_path(local_raw);
            if local_path.exists() {
                return Ok(canonical_or_display(&local_path));
            }
            if remote.is_none() {
                return Err(format!(
                    "Genome '{genome_id}' has a {kind}_local path '{}', but that file does not exist",
                    local_path.display()
                ));
            }
        }
        if let Some(remote_src) = remote {
            return Ok(remote_src.clone());
        }
        Err(format!(
            "Genome '{genome_id}' does not provide {kind}_local or {kind}_remote source"
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

fn infer_annotation_extension(source: &str) -> &'static str {
    let lower = source.to_ascii_lowercase();
    if lower.contains(".gff3") {
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

    let SourceReader { reader, total_bytes } = open_source_reader(source)?;
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
        let mut progress_reader = ProgressReader::new(reader, |done| on_progress(done, total_bytes));
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

fn parse_annotation_gene_records(path: &Path) -> Result<Vec<GenomeGeneRecord>, String> {
    let file = File::open(path)
        .map_err(|e| format!("Could not open annotation file '{}': {e}", path.display()))?;
    let reader = BufReader::new(file);
    let mut genes: Vec<GenomeGeneRecord> = vec![];
    for (line_no, line) in reader.lines().enumerate() {
        let line = line
            .map_err(|e| format!("Could not read annotation file '{}': {e}", path.display()))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        if !cols[2].eq_ignore_ascii_case("gene") {
            continue;
        }
        let start_1based = cols[3].parse::<usize>().map_err(|e| {
            format!(
                "Invalid gene start '{}' at line {} in '{}': {e}",
                cols[3],
                line_no + 1,
                path.display()
            )
        })?;
        let end_1based = cols[4].parse::<usize>().map_err(|e| {
            format!(
                "Invalid gene end '{}' at line {} in '{}': {e}",
                cols[4],
                line_no + 1,
                path.display()
            )
        })?;
        if start_1based == 0 || end_1based < start_1based {
            return Err(format!(
                "Invalid gene interval {}-{} at line {} in '{}'",
                start_1based,
                end_1based,
                line_no + 1,
                path.display()
            ));
        }
        let attrs = parse_annotation_attributes(cols[8]);
        let gene_name = attrs
            .get("gene_name")
            .or_else(|| attrs.get("name"))
            .or_else(|| attrs.get("gene"))
            .cloned();
        let gene_id = attrs
            .get("gene_id")
            .or_else(|| attrs.get("id"))
            .or_else(|| attrs.get("gene"))
            .cloned();
        let strand = cols[6].chars().next().and_then(|c| match c {
            '+' | '-' => Some(c),
            _ => None,
        });
        genes.push(GenomeGeneRecord {
            chromosome: cols[0].to_string(),
            start_1based,
            end_1based,
            strand,
            gene_id,
            gene_name,
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
    Ok(genes)
}

fn parse_annotation_attributes(raw: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for part in raw.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if let Some((key, value)) = part.split_once('=') {
            let key = key.trim().to_ascii_lowercase();
            let value = value.trim().trim_matches('"').to_string();
            if !key.is_empty() && !value.is_empty() {
                map.insert(key, value);
            }
            continue;
        }
        let mut pieces = part.splitn(2, char::is_whitespace);
        let key = pieces.next().unwrap_or_default().trim().to_ascii_lowercase();
        let value = pieces
            .next()
            .unwrap_or_default()
            .trim()
            .trim_matches('"')
            .trim()
            .to_string();
        if !key.is_empty() && !value.is_empty() {
            map.insert(key, value);
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
            "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
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
}
