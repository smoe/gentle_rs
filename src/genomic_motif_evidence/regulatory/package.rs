//! Strict, bounded metadata and selected-file verification for the regulatory
//! subset. No atlas payloads, package SQL, or package-wide glob is evaluated.

use super::*;
use ring::digest::{Context, SHA256};
use serde::Deserialize;

const MAX_INVENTORY_BYTES: u64 = 64 * 1024 * 1024;
const MAX_INVENTORY_ROWS: usize = 250_000;
const MAX_METADATA_FILE_BYTES: u64 = 512 * 1024 * 1024;

#[derive(Debug, Deserialize)]
pub(super) struct Manifest {
    pub schema_version: u32,
    pub kind: String,
    pub state: String,
    pub genome_id: String,
    pub assembly: String,
    pub annotation_release: String,
    pub regulatory_release: String,
    pub promoter_definition_id: String,
    pub production_plan_sha256: String,
    pub source_commit: String,
    pub database: String,
    pub coordinate_mode: String,
    pub scope: String,
    pub score_selection: String,
    pub upstream_bp: u64,
    pub downstream_bp: u64,
    pub chromosomes: Vec<String>,
    pub motif_count: usize,
    pub requires_tp73: bool,
    pub rows: u64,
    pub parquet_bytes: u64,
    pub complete_for_declared_annotation_and_source_floors: bool,
    pub complete_genome_scan: bool,
    pub files: Vec<RegulatoryMotifFileBinding>,
}

#[derive(Debug, Deserialize)]
pub(super) struct Annotation {
    pub schema_version: u32,
    pub kind: String,
    pub state: String,
    pub genome_id: String,
    pub assembly: String,
    pub annotation_release: String,
    pub regulatory_release: String,
    pub promoter_definition_id: String,
    pub upstream_bp: u64,
    pub downstream_bp: u64,
    pub coordinate_rule: String,
    pub regulatory_gff_sha256: String,
    pub coordinate_audit: Value,
    pub chromosomes: Vec<Chromosome>,
    pub files: Vec<RegulatoryMotifFileBinding>,
}

#[derive(Debug, Deserialize)]
pub(super) struct Chromosome {
    pub chrom: String,
    pub length: u64,
    pub feature_count: u64,
    pub promoter_count: u64,
    pub coverage: String,
    pub features: RegulatoryMotifFileBinding,
    pub promoters: RegulatoryMotifFileBinding,
}

#[derive(Debug, Clone, Deserialize)]
pub(super) struct Entry {
    pub chrom: String,
    pub motif_id: String,
    pub rows: u64,
    pub bytes: u64,
    pub path: Option<String>,
    pub sha256: Option<String>,
    pub state: String,
}

pub(super) struct Package {
    pub manifest: Manifest,
    pub annotation: Annotation,
    pub entries: BTreeMap<(String, String), Entry>,
    pub verified: Vec<RegulatoryMotifFileBinding>,
}

fn hash_valid(hash: &str) -> bool {
    hash.len() == 64
        && hash
            .bytes()
            .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase())
}

fn relative_path(value: &str) -> Result<(), String> {
    if value.is_empty()
        || value.contains(['\\', ':'])
        || value
            .split('/')
            .any(|p| p.is_empty() || p == "." || p == "..")
        || Path::new(value)
            .components()
            .any(|p| !matches!(p, std::path::Component::Normal(_)))
    {
        return Err(format!("Unsafe package-relative path: {value}"));
    }
    Ok(())
}

pub(super) fn checked_file(root: &Path, relative: &str) -> Result<PathBuf, String> {
    relative_path(relative)?;
    let mut path = root.to_path_buf();
    for part in relative.split('/') {
        path.push(part);
        if fs::symlink_metadata(&path)
            .map_err(|e| format!("{}: {e}", path.display()))?
            .file_type()
            .is_symlink()
        {
            return Err(format!("Symlink in package path: {relative}"));
        }
    }
    if !path.is_file() {
        return Err(format!("Not a package file: {relative}"));
    }
    Ok(path)
}

fn validate_bindings(files: &[RegulatoryMotifFileBinding]) -> Result<(), String> {
    let mut names = BTreeSet::new();
    for file in files {
        relative_path(&file.path)?;
        if !hash_valid(&file.sha256) || !names.insert(&file.path) {
            return Err("Malformed digest or duplicate file binding".into());
        }
    }
    Ok(())
}

pub(super) fn binding<'a>(
    files: &'a [RegulatoryMotifFileBinding],
    path: &str,
) -> Result<&'a RegulatoryMotifFileBinding, String> {
    files
        .iter()
        .find(|f| f.path == path)
        .ok_or_else(|| format!("Missing file binding: {path}"))
}

pub(super) fn verify(
    root: &Path,
    record: &RegulatoryMotifFileBinding,
    runtime: &Runtime,
    limit: u64,
) -> Result<PathBuf, String> {
    let path = checked_file(root, &record.path)?;
    if !hash_valid(&record.sha256) || record.bytes > limit {
        return Err(format!(
            "Invalid digest or file exceeds {limit} byte bound: {}",
            record.path
        ));
    }
    let mut file = fs::File::open(&path).map_err(|e| e.to_string())?;
    if file.metadata().map_err(|e| e.to_string())?.len() != record.bytes {
        return Err(format!("File size mismatch: {}", record.path));
    }
    let mut hash = Context::new(&SHA256);
    let mut buffer = [0; 64 * 1024];
    let mut read = 0_u64;
    loop {
        runtime.remaining()?;
        let count = file.read(&mut buffer).map_err(|e| e.to_string())?;
        if count == 0 {
            break;
        }
        read += count as u64;
        if read > record.bytes {
            return Err("File grew during verification".into());
        }
        hash.update(&buffer[..count]);
    }
    let digest = hash
        .finish()
        .as_ref()
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect::<String>();
    if read != record.bytes || digest != record.sha256 {
        return Err(format!("SHA-256 mismatch: {}", record.path));
    }
    Ok(path)
}

pub(super) fn json<T: serde::de::DeserializeOwned>(path: &Path, limit: u64) -> Result<T, String> {
    let file = fs::File::open(path).map_err(|e| e.to_string())?;
    if file.metadata().map_err(|e| e.to_string())?.len() > limit {
        return Err("JSON input exceeds size bound".into());
    }
    let mut bytes = Vec::new();
    file.take(limit + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| e.to_string())?;
    if bytes.len() as u64 > limit {
        return Err("JSON input grew beyond size bound".into());
    }
    serde_json::from_slice(&bytes).map_err(|e| e.to_string())
}

impl Package {
    pub fn open(paths: &PackagePaths, runtime: &Runtime) -> Result<Self, String> {
        checked_file(&paths.root, "manifest.json")?;
        let m: Manifest =
            serde_json::from_value(paths.manifest.clone()).map_err(|e| e.to_string())?;
        if m.schema_version != 1
            || m.kind != REGULATORY_MOTIF_PROVIDER
            || m.state != "complete"
            || m.coordinate_mode != "bed_0based_half_open"
            || m.scope != "regulatory_and_tss"
            || m.score_selection != "source_retention"
            || m.requires_tp73
            || m.complete_genome_scan
            || !m.complete_for_declared_annotation_and_source_floors
            || !hash_valid(&m.production_plan_sha256)
            || [
                &m.genome_id,
                &m.assembly,
                &m.annotation_release,
                &m.regulatory_release,
                &m.source_commit,
            ]
            .iter()
            .any(|s| s.is_empty())
            || m.upstream_bp > 10_000_000
            || m.downstream_bp > 10_000_000
            || m.motif_count == 0
            || m.chromosomes.is_empty()
            || m.chromosomes.iter().collect::<BTreeSet<_>>().len() != m.chromosomes.len()
            || m.motif_count
                .checked_mul(m.chromosomes.len())
                .is_none_or(|n| n > MAX_INVENTORY_ROWS)
        {
            return Err("Unsupported or inconsistent regulatory subset manifest".into());
        }
        if m.promoter_definition_id
            != format!(
                "tss_upstream_{}_downstream_{}_v1",
                m.upstream_bp, m.downstream_bp
            )
        {
            return Err("Promoter definition disagrees with declared flanks".into());
        }
        validate_bindings(&m.files)?;
        let mut verified = vec![];
        for name in [
            m.database.as_str(),
            "file_inventory.json",
            "annotation/manifest.json",
            "schema.sql",
        ] {
            let record = binding(&m.files, name)?;
            let path = verify(&paths.root, record, runtime, MAX_METADATA_FILE_BYTES)?;
            if name == m.database
                && path.canonicalize().map_err(|e| e.to_string())? != paths.database_path
            {
                return Err("Database override does not match the manifest-bound catalog".into());
            }
            verified.push(record.clone());
        }
        let annotation: Annotation = json(
            &paths.root.join("annotation/manifest.json"),
            MAX_MANIFEST_BYTES,
        )?;
        if annotation.schema_version != 1
            || annotation.kind != "regulatory_tfbs_annotation"
            || annotation.state != "complete"
            || annotation.genome_id != m.genome_id
            || annotation.assembly != m.assembly
            || annotation.annotation_release != m.annotation_release
            || annotation.regulatory_release != m.regulatory_release
            || annotation.promoter_definition_id != m.promoter_definition_id
            || annotation.upstream_bp != m.upstream_bp
            || annotation.downstream_bp != m.downstream_bp
            || annotation.coordinate_rule
                != "BED_half_open_offsets_include_TSS_base_clamped_to_sequence"
            || !hash_valid(&annotation.regulatory_gff_sha256)
            || annotation.coordinate_audit["status"] != "passed"
            || annotation.coordinate_audit["assembly"] != m.assembly
            || annotation.coordinate_audit["gff_sha256"] != annotation.regulatory_gff_sha256
        {
            return Err("Annotation identity/audit disagrees with regulatory subset".into());
        }
        validate_bindings(&annotation.files)?;
        let mut chromosomes = BTreeMap::new();
        for c in &annotation.chromosomes {
            if c.length == 0
                || c.length > i64::MAX as u64
                || !m.chromosomes.contains(&c.chrom)
                || chromosomes.insert(c.chrom.clone(), c).is_some()
                || !matches!(
                    c.coverage.as_str(),
                    "annotated_intersection" | "known_empty_intersection"
                )
                || (c.feature_count > 0 && c.promoter_count == 0)
                || (c.coverage == "annotated_intersection"
                    && (c.feature_count == 0 || c.promoter_count == 0))
                || binding(&annotation.files, &c.features.path)? != &c.features
                || binding(&annotation.files, &c.promoters.path)? != &c.promoters
            {
                return Err("Invalid or unsupported chromosome annotation coverage".into());
            }
        }
        if chromosomes.len() != m.chromosomes.len() {
            return Err("Missing chromosome annotation coverage".into());
        }
        let rows: Vec<Entry> = json(&paths.root.join("file_inventory.json"), MAX_INVENTORY_BYTES)?;
        if rows.len() != m.motif_count * m.chromosomes.len() {
            return Err("Incomplete chromosome/motif inventory".into());
        }
        let mut entries = BTreeMap::new();
        let mut motifs = BTreeSet::new();
        let mut payloads = BTreeSet::new();
        let (mut count, mut bytes) = (0_u64, 0_u64);
        for row in rows {
            runtime.remaining()?;
            let c = chromosomes
                .get(&row.chrom)
                .ok_or("Inventory chromosome absent from annotation")?;
            if row.motif_id.is_empty() {
                return Err("Empty inventory motif".into());
            }
            if row.state == "known_empty_intersection" {
                if c.coverage != row.state
                    || row.rows != 0
                    || row.bytes != 0
                    || row.path.is_some()
                    || row.sha256.is_some()
                {
                    return Err("Inconsistent known-empty intersection".into());
                }
            } else if row.state == "complete" && c.coverage == "annotated_intersection" {
                let path = row.path.as_deref().ok_or("Missing payload path")?;
                relative_path(path)?;
                if row.bytes == 0
                    || !row.sha256.as_deref().is_some_and(hash_valid)
                    || !payloads.insert(path.to_string())
                {
                    return Err("Invalid/duplicate inventory payload".into());
                }
            } else {
                return Err("Invalid inventory completion state".into());
            }
            count = count
                .checked_add(row.rows)
                .ok_or("Inventory count overflow")?;
            bytes = bytes
                .checked_add(row.bytes)
                .ok_or("Inventory size overflow")?;
            motifs.insert(row.motif_id.clone());
            if entries
                .insert((row.chrom.clone(), row.motif_id.clone()), row)
                .is_some()
            {
                return Err("Duplicate chromosome/motif inventory entry".into());
            }
        }
        if motifs.len() != m.motif_count || count != m.rows || bytes != m.parquet_bytes {
            return Err("Inventory totals disagree with manifest".into());
        }
        for chrom in &m.chromosomes {
            for motif in &motifs {
                if !entries.contains_key(&(chrom.clone(), motif.clone())) {
                    return Err("Incomplete inventory cross product".into());
                }
            }
        }
        Ok(Self {
            manifest: m,
            annotation,
            entries,
            verified,
        })
    }

    pub fn annotation_file(
        &mut self,
        paths: &PackagePaths,
        name: &str,
        runtime: &Runtime,
    ) -> Result<PathBuf, String> {
        let mut file = binding(&self.annotation.files, name)?.clone();
        file.path = format!("annotation/{}", file.path);
        let path = verify(&paths.root, &file, runtime, MAX_METADATA_FILE_BYTES)?;
        if !self.verified.contains(&file) {
            self.verified.push(file);
        }
        Ok(path)
    }
}
