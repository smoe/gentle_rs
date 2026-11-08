//! Bounded, explicitly configured local samtools selection; not an aligner or BAM decoder.
//! Full alignment hashes remain declarations unless an explicit full pass was requested.

use super::*;
use crate::engine::{PrimerPayloadFile, ReferenceContigDictionary};
use crate::supervised_process;
use std::{
    fs::{self, File},
    io::{BufRead, BufReader, Read, Write},
    path::{Path, PathBuf},
    process::Command,
    time::Duration,
};

pub const SCHEMA: &str = "gentle.indexed_locus_request.v1";
pub const RECEIPT_SCHEMA: &str = "gentle.indexed_locus_receipt.v1";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Format {
    Bam,
    Cram,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub schema: String,
    pub samtools: Option<String>,
    pub format: Format,
    pub alignment: PrimerPayloadFile,
    pub index: PrimerPayloadFile,
    /// Explicit producer declaration; a region pass cannot prove this pairing by full traversal.
    pub index_for_alignment_sha256: String,
    pub dictionary: PrimerPayloadFile,
    pub reference: ReferenceBinding,
    pub reference_id: String,
    pub reference_fasta: Option<PrimerPayloadFile>,
    pub reference_fai: Option<PrimerPayloadFile>,
    pub contig: String,
    pub interval: Interval,
    pub minimum_mapq: Option<u8>,
    pub require_flags: Option<u16>,
    pub exclude_flags: Option<u16>,
    pub verify_full_alignment: bool,
    pub timeout_ms: u64,
    pub max_sam_bytes: usize,
    pub max_stderr_bytes: usize,
    pub max_units: usize,
    pub text_limits: sam_text::Limits,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Receipt {
    pub schema: String,
    pub executable: PrimerPayloadFile,
    pub versions: String,
    pub argv: Vec<String>,
    pub environment: BTreeMap<String, String>,
    pub format: Format,
    pub alignment: PrimerPayloadFile,
    pub alignment_hash_verified: bool,
    pub index: PrimerPayloadFile,
    pub index_for_alignment_sha256_declaration: String,
    pub dictionary: PrimerPayloadFile,
    pub reference_fasta: Option<PrimerPayloadFile>,
    pub reference_fai: Option<PrimerPayloadFile>,
    pub contig: String,
    pub interval: Interval,
    pub timeout_ms: u64,
    pub max_sam_bytes: usize,
    pub max_stderr_bytes: usize,
    /// Hash of the emitted SAM header bytes, not the binary BAM header encoding.
    pub emitted_header_sha256: String,
    pub sam_sha256: String,
    pub sam_bytes: u64,
    pub stderr_bytes: u64,
    pub search_complete: bool,
    pub views_complete: bool,
    pub identity_scope: String,
}

impl Receipt {
    pub(crate) fn validate(&self, sam: &PrimerPayloadFile) -> Result<(), EngineError> {
        if self.schema != RECEIPT_SCHEMA
            || self.search_complete
            || self.views_complete
            || self.sam_sha256 != sam.sha256
            || self.sam_bytes != sam.bytes
            || self.index_for_alignment_sha256_declaration != self.alignment.sha256
            || self.identity_scope != "sam_declared_query_unit"
            || self.argv.is_empty()
            || !self.versions.contains("samtools")
            || !self.versions.contains("htslib")
            || self.interval.0 >= self.interval.1
        {
            return Err(fail("receipt"));
        }
        Ok(())
    }
}

pub struct Selected {
    _directory: tempfile::TempDir,
    pub sam: PrimerPayloadFile,
    pub manifest: PrimerPayloadFile,
    pub receipt: Receipt,
}

fn fail(reason: &str) -> EngineError {
    EngineError::invalid_input(format!("indexed_locus_{reason}"))
}
fn hash_ok(value: &str) -> bool {
    value.strip_prefix("sha256:").is_some_and(|v| {
        v.len() == 64
            && v.bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    })
}
fn io(e: impl std::fmt::Display) -> EngineError {
    fail(&e.to_string())
}

/// Validate raw spelling before any normalization and reject symlink/reparse-file inputs.
pub(crate) fn local_file(raw: &str) -> Result<PathBuf, EngineError> {
    let path = sam_intake::path(raw)?;
    if !path.is_absolute() || raw.contains("##idx##") || raw.contains("://") {
        return Err(fail("absolute_local_path"));
    }
    for parent in path.ancestors() {
        let meta = fs::symlink_metadata(parent).map_err(io)?;
        if meta.file_type().is_symlink() {
            return Err(fail("symlink"));
        }
        #[cfg(windows)]
        {
            use std::os::windows::fs::MetadataExt;
            if meta.file_attributes() & 0x400 != 0 {
                return Err(fail("reparse_point"));
            }
        }
    }
    if !fs::metadata(path).map_err(io)?.is_file() {
        return Err(fail("regular_file"));
    }
    Ok(path.to_owned())
}

struct Identity {
    owner: same_file::Handle,
    size: u64,
    modified: Option<std::time::SystemTime>,
}
fn identity(file: &PrimerPayloadFile) -> Result<Identity, EngineError> {
    let p = local_file(&file.path)?;
    let m = fs::metadata(&p).map_err(io)?;
    if !hash_ok(&file.sha256) || m.len() != file.bytes {
        return Err(fail("file_binding"));
    }
    Ok(Identity {
        owner: same_file::Handle::from_path(p).map_err(io)?,
        size: m.len(),
        modified: m.modified().ok(),
    })
}
fn unchanged(file: &PrimerPayloadFile, before: &Identity) -> Result<(), EngineError> {
    let after = identity(file)?;
    if before.owner != after.owner || before.size != after.size || before.modified != after.modified
    {
        return Err(fail("input_changed"));
    }
    Ok(())
}

pub(crate) fn fingerprint(
    path: &Path,
    cap: u64,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<PrimerPayloadFile, EngineError> {
    let raw = path.to_str().ok_or_else(|| fail("path_encoding"))?;
    let p = local_file(raw)?;
    let meta = fs::metadata(&p).map_err(io)?;
    if meta.len() > cap {
        return Err(fail("hash_byte_budget"));
    }
    let mut source = File::open(&p).map_err(io)?;
    let owner = same_file::Handle::from_file(source.try_clone().map_err(io)?).map_err(io)?;
    let mut context = ring::digest::Context::new(&ring::digest::SHA256);
    let (mut count, mut buffer) = (0u64, [0u8; 65536]);
    loop {
        check()?;
        let n = source.read(&mut buffer).map_err(io)?;
        if n == 0 {
            break;
        }
        count = count
            .checked_add(n as u64)
            .filter(|n| *n <= cap && *n <= meta.len())
            .ok_or_else(|| fail("hash_byte_budget"))?;
        context.update(&buffer[..n]);
    }
    let after = fs::metadata(&p).map_err(io)?;
    if count != meta.len()
        || after.len() != count
        || after.modified().ok() != meta.modified().ok()
        || same_file::Handle::from_path(&p).map_err(io)? != owner
    {
        return Err(fail("input_changed"));
    }
    Ok(PrimerPayloadFile {
        path: raw.into(),
        bytes: count,
        sha256: format!(
            "sha256:{}",
            context
                .finish()
                .as_ref()
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect::<String>()
        ),
    })
}
pub(crate) fn verify(
    file: &PrimerPayloadFile,
    cap: u64,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<(), EngineError> {
    let actual = fingerprint(Path::new(&file.path), cap, check)?;
    if actual.bytes != file.bytes || actual.sha256 != file.sha256 {
        return Err(fail("hash_binding"));
    }
    Ok(())
}
fn write_json<T: Serialize>(
    path: &Path,
    value: &T,
    cap: u64,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<PrimerPayloadFile, EngineError> {
    let mut f = File::options()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(io)?;
    serde_json::to_writer(&mut f, value).map_err(io)?;
    f.flush().map_err(io)?;
    f.sync_all().map_err(io)?;
    drop(f);
    fingerprint(path, cap, check)
}

fn manifest(
    sam: &PrimerPayloadFile,
    dictionary: &ReferenceContigDictionary,
    request: &Request,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<(sam_intake::Manifest, String), EngineError> {
    let contigs: BTreeMap<_, _> = dictionary
        .contigs
        .iter()
        .map(|c| (c.name.clone(), c.length as usize))
        .collect();
    let mut decoder = sam_text::Consumer::with_version_policy(
        &contigs,
        &request.text_limits,
        sam_text::VersionPolicy::Strict16,
        check,
    )?;
    let mut source = BufReader::new(File::open(&sam.path).map_err(io)?);
    let mut header = ring::digest::Context::new(&ring::digest::SHA256);
    let mut in_header = true;
    loop {
        check()?;
        let mut line = vec![];
        let n = source
            .by_ref()
            .take(request.text_limits.max_line_bytes as u64 + 1)
            .read_until(b'\n', &mut line)
            .map_err(io)?;
        if n == 0 {
            break;
        }
        if in_header && line.starts_with(b"@") {
            header.update(&line);
        } else {
            in_header = false;
        }
        decoder.push_line(&line, check)?;
    }
    let decoded = decoder.finish(check)?;
    if decoded.raw_sha256 != sam.sha256 || decoded.raw_bytes as u64 != sam.bytes {
        return Err(fail("sam_binding"));
    }
    let mut grouped: BTreeMap<(String, Option<String>), Vec<sam_text::DecodedRow>> =
        BTreeMap::new();
    for row in decoded.rows {
        check()?;
        grouped
            .entry((row.row.qname.clone(), row.read_group.clone()))
            .or_default()
            .push(row);
    }
    if grouped.len() > request.max_units {
        return Err(fail("unit_budget"));
    }
    let mut units = vec![];
    for (index, (_, rows)) in grouped.into_iter().enumerate() {
        let axes: BTreeSet<_> = rows
            .iter()
            .map(|r| sam::normalized_axis_length(&r.row))
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flatten()
            .collect();
        if axes.len() > 1 {
            return Err(fail("query_axis_conflict"));
        }
        units.push(sam_intake::Unit {unit_id:format!("region_unit_{index:08}"),unclipped_axis_length:axes.into_iter().next(),
            rows:rows.into_iter().map(|r| sam_intake::RowWitness {ordinal:r.row.ordinal,qname:r.row.qname,read_group:r.read_group}).collect(),
            masks:None,ends:[vec![],vec![]],search:Some(sam_intake::SearchDeclaration {search_complete:false,views_complete:false,
                limitations:vec!["Indexed region selection; other alignments of each QNAME may be unseen. QNAME+RG grouping is per staged SAM, not original-read identity.".into()],
                producer:"GENtle configured samtools".into(),version:env!("CARGO_PKG_VERSION").into(),arguments:vec![],selection_scope:format!("{}:{}..{} half-open",request.contig,request.interval.0,request.interval.1)})});
    }
    Ok((
        sam_intake::Manifest {
            schema: "gentle.observed_sam_query_units.v1".into(),
            sam_sha256: sam.sha256.clone(),
            units,
        },
        format!(
            "sha256:{}",
            header
                .finish()
                .as_ref()
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect::<String>()
        ),
    ))
}

/// Returns owned, closed private staging files for the unchanged S2.3 intake. No publication.
pub fn select(
    request: &Request,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Selected, EngineError> {
    check()?;
    let executable = request
        .samtools
        .as_ref()
        .ok_or_else(|| fail("unavailable_samtools_not_configured"))?;
    if request.schema != SCHEMA
        || request.interval.0 >= request.interval.1
        || request.interval.1 > i32::MAX as usize
        || request.contig.is_empty()
        || request.contig.chars().any(char::is_whitespace)
        || !(1..=3_600_000).contains(&request.timeout_ms)
        || !(1..=128 * 1024 * 1024).contains(&request.max_sam_bytes)
        || !(1..=1024 * 1024).contains(&request.max_stderr_bytes)
        || !(1..=4096).contains(&request.max_units)
        || request.index_for_alignment_sha256 != request.alignment.sha256
        || request.dictionary.bytes > 16 * 1024 * 1024
        || request.index.bytes > 512 * 1024 * 1024
        || request
            .reference_fai
            .as_ref()
            .is_some_and(|f| f.bytes > 64 * 1024 * 1024)
        || request.text_limits.max_input_bytes > 128 * 1024 * 1024
        || request.text_limits.max_retained_bytes > 128 * 1024 * 1024
    {
        return Err(fail("admission"));
    }
    let executable = fingerprint(&local_file(executable)?, 128 * 1024 * 1024, check)?;
    let mut bindings = vec![
        request.alignment.clone(),
        request.index.clone(),
        request.dictionary.clone(),
    ];
    match request.format {
        Format::Cram => {
            let fasta = request
                .reference_fasta
                .as_ref()
                .ok_or_else(|| fail("cram_reference_required"))?;
            let fai = request
                .reference_fai
                .as_ref()
                .ok_or_else(|| fail("cram_fai_required"))?;
            if fasta.sha256 != request.reference.fasta_sha256
                || fai.path != format!("{}.fai", fasta.path)
                || fasta.path.ends_with(".gz")
            {
                return Err(fail("cram_reference_binding"));
            }
            bindings.extend([fasta.clone(), fai.clone()]);
        }
        Format::Bam if request.reference_fasta.is_some() || request.reference_fai.is_some() => {
            return Err(fail("unexpected_reference"));
        }
        _ => {}
    }
    let identities = bindings
        .iter()
        .map(identity)
        .collect::<Result<Vec<_>, _>>()?;
    for file in &bindings[1..] {
        verify(file, file.bytes, check)?;
    }
    if request.verify_full_alignment {
        verify(&request.alignment, request.alignment.bytes, check)?;
    }
    if identities[1]
        .modified
        .zip(identities[0].modified)
        .is_none_or(|(index, data)| index < data)
    {
        return Err(fail("stale_index"));
    }
    let dictionary: ReferenceContigDictionary =
        serde_json::from_reader(File::open(&request.dictionary.path).map_err(io)?).map_err(io)?;
    if request.dictionary.bytes > 16 * 1024 * 1024
        || dictionary.reference_id != request.reference_id
        || dictionary.assembly != request.reference.assembly
        || dictionary.reference_fasta_sha256 != request.reference.fasta_sha256
        || !dictionary
            .contigs
            .iter()
            .any(|c| c.name == request.contig && c.length >= request.interval.1 as u64)
    {
        return Err(fail("dictionary"));
    }
    if let Some(fai) = &request.reference_fai {
        let lengths = crate::genomes::fasta_index_lengths(Path::new(&fai.path)).map_err(io)?;
        if dictionary
            .contigs
            .iter()
            .any(|c| lengths.get(&c.name) != Some(&c.length))
        {
            return Err(fail("cram_fai_dictionary"));
        }
    }
    let directory = tempfile::tempdir().map_err(io)?;
    // Canonical private directory avoids platform aliases while raw caller paths stay checked first.
    let root = directory.path().canonicalize().map_err(io)?;
    let environment = BTreeMap::from([
        (
            "REF_PATH".into(),
            root.join("no-reference-lookup/%s")
                .to_string_lossy()
                .into_owned(),
        ),
        (
            "REF_CACHE".into(),
            root.join("no-reference-cache/%s")
                .to_string_lossy()
                .into_owned(),
        ),
        ("LC_ALL".into(), "C".into()),
    ]);
    fs::create_dir(root.join("no-reference-lookup")).map_err(io)?;
    fs::create_dir(root.join("no-reference-cache")).map_err(io)?;
    let command = |args: &[String]| {
        let mut c = Command::new(&executable.path);
        c.args(args).envs(&environment);
        c
    };
    let limits = supervised_process::Limits {
        timeout: Duration::from_millis(request.timeout_ms),
        stdout_bytes: request.max_sam_bytes,
        stderr_bytes: request.max_stderr_bytes,
        refuse_overflow: true,
    };
    let mut cancelled = false;
    let version = supervised_process::run(
        &mut command(&["--version".into()]),
        "samtools version",
        supervised_process::Limits {
            stdout_bytes: 65536,
            ..limits
        },
        None,
        &mut || {
            cancelled |= check().is_err();
            cancelled
        },
    )
    .map_err(io)?;
    let versions = String::from_utf8(version.output.stdout).map_err(io)?;
    if !version.output.status.success()
        || !versions.contains("samtools")
        || !versions.contains("htslib")
    {
        return Err(fail("tool_version"));
    }
    // Read only the header first. Reject CRAM UR fallback locations before decoding any slice.
    let header_args = vec![
        "view".into(),
        "-H".into(),
        "--no-PG".into(),
        request.alignment.path.clone(),
    ];
    let header = supervised_process::run(
        &mut command(&header_args),
        "samtools header",
        supervised_process::Limits {
            stdout_bytes: 1024 * 1024,
            ..limits
        },
        None,
        &mut || {
            cancelled |= check().is_err();
            cancelled
        },
    )
    .map_err(io)?;
    if !header.output.status.success() || !header.output.stderr.is_empty() {
        return Err(fail("header_refused"));
    }
    if request.format == Format::Cram {
        let header = std::str::from_utf8(&header.output.stdout).map_err(io)?;
        let reference = request
            .reference_fasta
            .as_ref()
            .ok_or_else(|| fail("cram_reference_required"))?;
        for uri in header
            .lines()
            .filter(|line| line.starts_with("@SQ\t"))
            .flat_map(|line| line.split('\t'))
            .filter_map(|field| field.strip_prefix("UR:"))
        {
            if uri != reference.path {
                return Err(fail("cram_reference_fallback_refused"));
            }
        }
    }
    let bed = root.join("region.bed");
    {
        let mut file = File::options()
            .write(true)
            .create_new(true)
            .open(&bed)
            .map_err(io)?;
        writeln!(
            file,
            "{}\t{}\t{}",
            request.contig, request.interval.0, request.interval.1
        )
        .map_err(io)?;
        file.sync_all().map_err(io)?;
    }
    let mut argv = vec![
        "view".into(),
        "-h".into(),
        "--no-PG".into(),
        "-M".into(),
        "-L".into(),
        bed.to_string_lossy().into_owned(),
        "-X".into(),
    ];
    for (flag, value) in [
        ("-q", request.minimum_mapq.map(u32::from)),
        ("-f", request.require_flags.map(u32::from)),
        ("-F", request.exclude_flags.map(u32::from)),
    ] {
        if let Some(value) = value {
            argv.extend([flag.into(), value.to_string()]);
        }
    }
    if let Some(reference) = &request.reference_fasta {
        argv.extend(["-T".into(), reference.path.clone()]);
    }
    argv.extend([request.alignment.path.clone(), request.index.path.clone()]);
    let sam_path = root.join("selection.sam");
    let file = File::options()
        .write(true)
        .create_new(true)
        .open(&sam_path)
        .map_err(io)?;
    let result = supervised_process::run(
        &mut command(&argv),
        "samtools selection",
        limits,
        Some(file),
        &mut || {
            cancelled |= check().is_err();
            cancelled
        },
    )
    .map_err(io)?;
    check()?;
    // Warnings include stale/wrong index and reference problems: refuse all stderr conservatively.
    if !result.output.status.success() || !result.output.stderr.is_empty() {
        return Err(fail(&format!(
            "tool_refused: {}",
            String::from_utf8_lossy(&result.output.stderr)
        )));
    }
    for (file, before) in bindings.iter().zip(&identities) {
        unchanged(file, before)?;
    }
    verify(&executable, executable.bytes, check)?;
    let sam = fingerprint(&sam_path, request.max_sam_bytes as u64, check)?;
    if sam.bytes != result.stdout_bytes {
        return Err(fail("staged_output_binding"));
    }
    let (manifest, header) = manifest(&sam, &dictionary, request, check)?;
    let manifest = write_json(
        &root.join("query-units.json"),
        &manifest,
        128 * 1024 * 1024,
        check,
    )?;
    let receipt = Receipt {
        schema: RECEIPT_SCHEMA.into(),
        executable,
        versions,
        argv,
        environment,
        format: request.format,
        alignment: request.alignment.clone(),
        alignment_hash_verified: request.verify_full_alignment,
        index: request.index.clone(),
        index_for_alignment_sha256_declaration: request.index_for_alignment_sha256.clone(),
        dictionary: request.dictionary.clone(),
        reference_fasta: request.reference_fasta.clone(),
        reference_fai: request.reference_fai.clone(),
        contig: request.contig.clone(),
        interval: request.interval,
        timeout_ms: request.timeout_ms,
        max_sam_bytes: request.max_sam_bytes,
        max_stderr_bytes: request.max_stderr_bytes,
        emitted_header_sha256: header,
        sam_sha256: sam.sha256.clone(),
        sam_bytes: sam.bytes,
        stderr_bytes: result.stderr_bytes,
        search_complete: false,
        views_complete: false,
        identity_scope: "sam_declared_query_unit".into(),
    };
    receipt.validate(&sam)?;
    Ok(Selected {
        _directory: directory,
        sam,
        manifest,
        receipt,
    })
}

#[cfg(test)]
#[path = "indexed_locus_tests.rs"]
mod tests;
