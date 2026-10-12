//! Bounded, read-only validation of versioned TSS FASTA bundles.
//!
//! The canonical `gentle.tss_fasta_bundle.v1` checksum list binds its manifest.
//! The explicit `gentle.target_tss_fasta_export.v1` adapter instead verifies the
//! manifest's checksum-list digest, then its FASTA members. Both hash original
//! file bytes and canonical uppercase A/C/G/T/N sequences separately. Sequences
//! are already transcript-oriented, including minus-strand records. Verification
//! establishes internal consistency, not independent reference authenticity.

use crate::digest_utils::sha256_hex_bytes;
use gentle_engine::tss_profiles::{
    validate_manifest, validate_record, validate_selection, validate_sha256,
};
use gentle_protocol::EngineError;
use gentle_protocol::tss_profiles::{
    BUNDLE_SCHEMA, ComputeTssProfilesRequest, TssBundleManifest, TssBundleSource, TssGeometry,
    TssInputBinding, TssRecord, TssReference, TssSelection, TssSelectionEvidence, TssStrand,
};
use serde::de::{self, DeserializeOwned, MapAccess, SeqAccess, Visitor};
use serde::{Deserialize, Deserializer};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::fs::{self, File};
use std::io::Read;
use std::path::{Path, PathBuf};

mod target_export;

const CHECKSUMS_NAME: &str = "SHA256SUMS";
const MAX_MANIFEST_BYTES: u64 = 8 * 1024 * 1024;
const MAX_CHECKSUM_BYTES: u64 = 1024 * 1024;
const MAX_MEMBER_BYTES: u64 = 32 * 1024 * 1024;
const MAX_TOTAL_BYTES: u64 = 64 * 1024 * 1024;
const MAX_CHECKSUM_ENTRIES: usize = 4_096;
const MAX_HEADER_BYTES: usize = 64 * 1024;
const MAX_SEQUENCE_BASES: usize = 1_000_000;
const ORIENTATION: &str = "transcript_5prime_to_3prime";
const HEADER_KEYS: [&str; 12] = [
    "gene_symbol",
    "gene_id",
    "promoter_id",
    "assembly",
    "chromosome",
    "strand",
    "tss_1based",
    "genomic_1based",
    "window",
    "orientation",
    "transcripts",
    "sequence_sha256",
];

/// Internally consistent input windows, with portable byte bindings.
#[derive(Debug)]
pub struct VerifiedTssBundle {
    pub reference: TssReference,
    /// Identity of the original manifest, not an adapter-generated representation.
    pub source: TssBundleSource,
    /// Manifest order, canonical transcript-oriented sequence, explicit selection.
    pub records: Vec<(TssRecord, String, bool)>,
    pub inputs: Vec<TssInputBinding>,
    pub warnings: Vec<String>,
    /// Descriptive source evidence keyed only by the joined promoter ID.
    pub selection_evidence: BTreeMap<String, TssSelectionEvidence>,
}

struct BundleFiles {
    root: PathBuf,
    declared_parent: PathBuf,
    manifest_name: String,
    manifest_path: PathBuf,
    manifest_bytes: Vec<u8>,
    sums_path: PathBuf,
    sums_bytes: Vec<u8>,
    sums: BTreeMap<String, String>,
    total_bytes: u64,
}

/// Read and validate all inputs before returning any usable bundle records.
///
/// An explicit FASTA list must name exactly the manifest's contained members.
/// Canonical selections are contained checksum-bound members. Target-export
/// selections may also be separately supplied absolute local files; their exact
/// bytes are receipt-bound, never discovered by following manifest links.
/// No panel resolution, downloads, reference extraction or output writes occur.
pub fn read_bundle(request: &ComputeTssProfilesRequest) -> Result<VerifiedTssBundle, EngineError> {
    let manifest_path = Path::new(&request.manifest);
    let manifest_name = manifest_path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| EngineError::invalid_input("Manifest needs a UTF-8 filename"))?;
    validate_relative_name(manifest_name)?;
    if manifest_name == CHECKSUMS_NAME {
        return Err(EngineError::invalid_input(
            "Manifest and SHA256SUMS must be distinct files",
        ));
    }
    let parent = manifest_path
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let root = fs::canonicalize(parent).map_err(|error| {
        EngineError::invalid_input("Cannot resolve the TSS bundle directory").with_cause(error)
    })?;
    let mut total_bytes = 0;
    let manifest_path = contained_file(&root, manifest_name)?;
    let manifest_bytes = read_bounded(&manifest_path, MAX_MANIFEST_BYTES, &mut total_bytes)?;
    let sums_path = contained_file(&root, CHECKSUMS_NAME)?;
    if sums_path == manifest_path {
        return Err(EngineError::invalid_input(
            "Manifest and SHA256SUMS cannot alias the same file",
        ));
    }
    let sums_bytes = read_bounded(&sums_path, MAX_CHECKSUM_BYTES, &mut total_bytes)?;
    let sums = parse_checksums(&sums_bytes)?;
    #[derive(Deserialize)]
    struct Schema {
        schema: String,
    }
    let schema: Schema = parse_json(&manifest_bytes, "bundle manifest")?;
    if schema.schema == target_export::SCHEMA {
        return target_export::read_target_bundle(
            request,
            BundleFiles {
                root,
                declared_parent: parent.to_path_buf(),
                manifest_name: manifest_name.into(),
                manifest_path,
                manifest_bytes,
                sums_path,
                sums_bytes,
                sums,
                total_bytes,
            },
        );
    }
    if schema.schema != BUNDLE_SCHEMA {
        return Err(EngineError::invalid_input(format!(
            "Unsupported TSS bundle schema {:?}",
            schema.schema
        )));
    }
    verify_checksum(&sums, manifest_name, &manifest_bytes)?;
    let manifest: TssBundleManifest = parse_json(&manifest_bytes, "bundle manifest")?;
    validate_manifest(&manifest)?;
    check_reference_expectations(request, &manifest.reference)?;
    if request.expected_dataset_id.is_some() {
        return Err(EngineError::invalid_input(
            "Canonical bundle does not declare a dataset_id to match the requested expectation",
        ));
    }

    for record in &manifest.records {
        if record
            .geometry
            .length()
            .is_none_or(|n| n > MAX_SEQUENCE_BASES)
        {
            return Err(EngineError::invalid_input(format!(
                "Promoter {} exceeds the {MAX_SEQUENCE_BASES}-base window limit",
                record.promoter_id
            )));
        }
    }
    for (name, digest) in &manifest.fasta_files {
        validate_relative_name(name)?;
        if name == manifest_name || name == CHECKSUMS_NAME {
            return Err(EngineError::invalid_input(
                "FASTA members cannot reuse the manifest or checksum filename",
            ));
        }
        if sums.get(name) != Some(digest) {
            return Err(EngineError::invalid_input(format!(
                "FASTA manifest digest does not match SHA256SUMS for {name}"
            )));
        }
    }
    check_requested_fastas(request, &root, parent, &manifest.fasta_files)?;
    let selection_name = request
        .selection
        .as_deref()
        .map(|name| input_name(&root, parent, name))
        .transpose()?;
    if let Some(name) = &selection_name {
        if name == manifest_name
            || name == CHECKSUMS_NAME
            || manifest.fasta_files.contains_key(name)
            || !sums.contains_key(name)
        {
            return Err(EngineError::invalid_input(
                "Selection must be a distinct bundle member listed in SHA256SUMS",
            ));
        }
    }

    let mut inputs = vec![
        binding("bundle_manifest", manifest_name, &manifest_bytes),
        binding("bundle_checksums", CHECKSUMS_NAME, &sums_bytes),
    ];
    let expected: BTreeMap<_, _> = manifest
        .records
        .iter()
        .map(|record| (record.promoter_id.as_str(), record))
        .collect();
    let mut sequences = BTreeMap::new();
    let mut selected = BTreeSet::new();
    let mut resolved_paths = BTreeSet::from([manifest_path, sums_path]);
    for (name, digest) in &sums {
        if name == manifest_name {
            continue;
        }
        let path = contained_file(&root, name)?;
        if !resolved_paths.insert(path.clone()) {
            return Err(EngineError::invalid_input(format!(
                "Bundle paths alias the same file: {name}"
            )));
        }
        let limit = if selection_name.as_ref() == Some(name) {
            MAX_MANIFEST_BYTES
        } else {
            MAX_MEMBER_BYTES
        };
        let bytes = read_bounded(&path, limit, &mut total_bytes)?;
        let actual_digest = sha256_hex_bytes(&bytes);
        if &actual_digest != digest {
            return Err(EngineError::invalid_input(format!(
                "SHA256SUMS mismatch for {name}"
            )));
        }
        let role = if manifest.fasta_files.contains_key(name) {
            parse_fasta(&bytes, name, &manifest.reference, &expected, &mut sequences)?;
            "fasta"
        } else if selection_name.as_ref() == Some(name) {
            let selection: TssSelection = parse_json(&bytes, "TSS selection")?;
            validate_selection(&selection, &manifest)?;
            selected.extend(
                selection
                    .selected
                    .into_iter()
                    .map(|record| record.promoter_id),
            );
            "selection"
        } else {
            "bundle_auxiliary"
        };
        inputs.push(TssInputBinding {
            role: role.into(),
            name: name.clone(),
            sha256: actual_digest,
        });
    }
    finish_bundle(
        manifest.records,
        sequences,
        &selected,
        VerifiedTssBundle {
            reference: manifest.reference,
            source: TssBundleSource {
                schema: manifest.schema,
                manifest_sha256: sha256_hex_bytes(&manifest_bytes),
                source_revision: None,
                dataset_id: None,
                producer_sha256: None,
            },
            records: Vec::new(),
            inputs,
            warnings: Vec::new(),
            selection_evidence: BTreeMap::new(),
        },
    )
}

fn finish_bundle(
    manifest_records: Vec<TssRecord>,
    mut sequences: BTreeMap<String, String>,
    selected: &BTreeSet<String>,
    mut bundle: VerifiedTssBundle,
) -> Result<VerifiedTssBundle, EngineError> {
    if sequences.len() != manifest_records.len() {
        let missing: Vec<_> = manifest_records
            .iter()
            .filter(|record| !sequences.contains_key(&record.promoter_id))
            .map(|record| record.promoter_id.as_str())
            .take(10)
            .collect();
        return Err(EngineError::invalid_input(format!(
            "FASTA/manifest promoter membership mismatch; missing {}",
            missing.join(", ")
        )));
    }
    bundle.warnings.insert(0,
        "Bundle hashes and internal consistency verified only; reference genome/annotation authenticity and genomic sequence extraction were not independently verified. No downloads were performed.".into());
    let mut sequence_loci = BTreeMap::new();
    for record in manifest_records {
        if let Some(first_id) = sequence_loci.get(&record.sequence_sha256) {
            bundle.warnings.push(format!(
                "Repeated sequence SHA-256 {} at distinct loci (promoters {first_id} and {}); both positional identities retained",
                record.sequence_sha256, record.promoter_id
            ));
        } else {
            sequence_loci.insert(record.sequence_sha256.clone(), record.promoter_id.clone());
        }
        let sequence = sequences.remove(&record.promoter_id).ok_or_else(|| {
            EngineError::invalid_input("FASTA/manifest promoter membership mismatch")
        })?;
        let is_selected = selected.contains(&record.promoter_id);
        bundle.records.push((record, sequence, is_selected));
    }
    Ok(bundle)
}

fn check_requested_fastas(
    request: &ComputeTssProfilesRequest,
    root: &Path,
    parent: &Path,
    files: &BTreeMap<String, String>,
) -> Result<(), EngineError> {
    if !request.fasta.is_empty() {
        let mut requested = BTreeSet::new();
        for name in &request.fasta {
            let name = input_name(root, parent, name)?;
            if !requested.insert(name) {
                return Err(EngineError::invalid_input(
                    "Requested FASTA list contains duplicate paths",
                ));
            }
        }
        if requested != files.keys().cloned().collect() {
            return Err(EngineError::invalid_input(
                "Requested FASTA list must exactly match the manifest's listed paths",
            ));
        }
    }
    Ok(())
}

fn binding(role: &str, name: &str, bytes: &[u8]) -> TssInputBinding {
    TssInputBinding {
        role: role.into(),
        name: name.into(),
        sha256: sha256_hex_bytes(bytes),
    }
}

fn check_reference_expectations(
    request: &ComputeTssProfilesRequest,
    reference: &TssReference,
) -> Result<(), EngineError> {
    if request.expected_genome_id != reference.genome_id
        || request
            .expected_assembly
            .as_ref()
            .is_some_and(|assembly| assembly != &reference.assembly)
        || request
            .expected_annotation_release
            .as_ref()
            .is_some_and(|release| reference.annotation_release.as_ref() != Some(release))
    {
        return Err(EngineError::invalid_input(
            "Bundle reference does not exactly match the requested genome/assembly/annotation expectations",
        ));
    }
    Ok(())
}

fn validate_relative_name(name: &str) -> Result<(), EngineError> {
    if name.is_empty()
        || name.len() > 4_096
        || name.trim() != name
        || name
            .chars()
            .any(|c| c.is_control() || matches!(c, '\\' | ':'))
        || name
            .split('/')
            .any(|part| part.is_empty() || matches!(part, "." | ".."))
    {
        return Err(EngineError::invalid_input(format!(
            "Bundle path must be a portable relative filename without traversal: {name:?}"
        )));
    }
    Ok(())
}

fn input_name(root: &Path, declared_parent: &Path, value: &str) -> Result<String, EngineError> {
    let path = Path::new(value);
    let name = if path.is_absolute() {
        path.strip_prefix(root)
            .or_else(|_| path.strip_prefix(declared_parent))
            .ok()
            .and_then(Path::to_str)
            .ok_or_else(|| {
                EngineError::invalid_input("Requested input must be inside the bundle")
            })?
    } else {
        value
    };
    validate_relative_name(name)?;
    Ok(name.into())
}

fn contained_file(root: &Path, name: &str) -> Result<PathBuf, EngineError> {
    validate_relative_name(name)?;
    let path = fs::canonicalize(root.join(name)).map_err(|error| {
        EngineError::invalid_input(format!("Cannot resolve bundle member {name}")).with_cause(error)
    })?;
    if !path.starts_with(root) || path == root {
        return Err(EngineError::invalid_input(format!(
            "Bundle member escapes its directory, including through symlinks: {name}"
        )));
    }
    let metadata = fs::metadata(&path).map_err(|error| {
        EngineError::invalid_input(format!("Cannot inspect bundle member {name}")).with_cause(error)
    })?;
    if !metadata.is_file() {
        return Err(EngineError::invalid_input(format!(
            "Bundle member must be a regular file: {name}"
        )));
    }
    Ok(path)
}

/// Open only regular inputs; Unix nonblocking open also closes the FIFO-swap race.
pub(crate) fn open_regular_input(path: &Path) -> Result<File, EngineError> {
    let inspect = |metadata: fs::Metadata| {
        if metadata.is_file() {
            Ok(())
        } else {
            Err(EngineError::invalid_input(
                "TSS input must be a regular file",
            ))
        }
    };
    inspect(fs::metadata(path).map_err(|error| {
        EngineError::invalid_input("Cannot inspect TSS input file").with_cause(error)
    })?)?;
    let mut options = fs::OpenOptions::new();
    options.read(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        options.custom_flags(libc::O_NONBLOCK);
    }
    let file = options.open(path).map_err(|error| {
        EngineError::invalid_input("Cannot open TSS input file").with_cause(error)
    })?;
    inspect(file.metadata().map_err(|error| {
        EngineError::invalid_input("Cannot inspect opened TSS input").with_cause(error)
    })?)?;
    Ok(file)
}

fn read_bounded(path: &Path, limit: u64, total: &mut u64) -> Result<Vec<u8>, EngineError> {
    let limit = limit.min(MAX_TOTAL_BYTES.saturating_sub(*total));
    let file = open_regular_input(path)?;
    let metadata = file.metadata().map_err(|error| {
        EngineError::invalid_input("Cannot inspect bundle file").with_cause(error)
    })?;
    if !metadata.is_file() || metadata.len() > limit {
        return Err(EngineError::invalid_input(
            "Bundle file exceeds its byte limit or total input budget, or is not a regular file",
        ));
    }
    let mut bytes = Vec::new();
    file.take(limit + 1)
        .read_to_end(&mut bytes)
        .map_err(|error| EngineError::invalid_input("Cannot read bundle file").with_cause(error))?;
    if bytes.len() as u64 > limit {
        return Err(EngineError::invalid_input(
            "Bundle file grew beyond its byte limit or total input budget during reading",
        ));
    }
    *total += bytes.len() as u64;
    Ok(bytes)
}

fn parse_checksums(bytes: &[u8]) -> Result<BTreeMap<String, String>, EngineError> {
    let text = std::str::from_utf8(bytes)
        .map_err(|_| EngineError::invalid_input("SHA256SUMS must be UTF-8 text"))?;
    let mut sums = BTreeMap::new();
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let digest = line
            .get(..64)
            .ok_or_else(|| EngineError::invalid_input("Malformed SHA256SUMS digest"))?;
        validate_sha256(digest)?;
        if !matches!(line.get(64..66), Some("  " | " *")) {
            return Err(EngineError::invalid_input(
                "SHA256SUMS rows require '<lowercase digest>  <path>' or '<lowercase digest> *<path>'",
            ));
        }
        let name = &line[66..];
        validate_relative_name(name)?;
        if name == CHECKSUMS_NAME {
            return Err(EngineError::invalid_input(
                "SHA256SUMS must not include a self-referential checksum",
            ));
        }
        if sums.insert(name.into(), digest.into()).is_some() {
            return Err(EngineError::invalid_input(format!(
                "Duplicate SHA256SUMS path {name}"
            )));
        }
        if sums.len() > MAX_CHECKSUM_ENTRIES {
            return Err(EngineError::invalid_input("Too many SHA256SUMS members"));
        }
    }
    Ok(sums)
}

fn verify_checksum(
    sums: &BTreeMap<String, String>,
    name: &str,
    bytes: &[u8],
) -> Result<(), EngineError> {
    if sums.get(name) != Some(&sha256_hex_bytes(bytes)) {
        return Err(EngineError::invalid_input(format!(
            "SHA256SUMS missing or mismatched checksum for {name}"
        )));
    }
    Ok(())
}

fn parse_json<T: DeserializeOwned>(bytes: &[u8], label: &str) -> Result<T, EngineError> {
    serde_json::from_slice::<UniqueJsonKeys>(bytes).map_err(|error| {
        EngineError::invalid_input(format!("Invalid {label} JSON")).with_cause(error)
    })?;
    serde_json::from_slice(bytes).map_err(|error| {
        EngineError::invalid_input(format!("Invalid {label} schema")).with_cause(error)
    })
}

// Serde structs reject duplicate fields; maps otherwise silently keep the last
// filename. Walk all objects first so conflicting digest keys cannot be hidden.
struct UniqueJsonKeys;

impl<'de> Deserialize<'de> for UniqueJsonKeys {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct UniqueVisitor;
        impl<'de> Visitor<'de> for UniqueVisitor {
            type Value = UniqueJsonKeys;

            fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str("JSON with unique object keys")
            }

            fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
                let mut keys = BTreeSet::new();
                while let Some(key) = map.next_key::<String>()? {
                    if !keys.insert(key) {
                        return Err(de::Error::custom("duplicate JSON object key"));
                    }
                    map.next_value::<UniqueJsonKeys>()?;
                }
                Ok(UniqueJsonKeys)
            }

            fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
                while seq.next_element::<UniqueJsonKeys>()?.is_some() {}
                Ok(UniqueJsonKeys)
            }

            fn visit_bool<E: de::Error>(self, _: bool) -> Result<Self::Value, E> {
                Ok(UniqueJsonKeys)
            }

            fn visit_i64<E: de::Error>(self, _: i64) -> Result<Self::Value, E> {
                Ok(UniqueJsonKeys)
            }

            fn visit_u64<E: de::Error>(self, _: u64) -> Result<Self::Value, E> {
                Ok(UniqueJsonKeys)
            }

            fn visit_f64<E: de::Error>(self, _: f64) -> Result<Self::Value, E> {
                Ok(UniqueJsonKeys)
            }

            fn visit_str<E: de::Error>(self, _: &str) -> Result<Self::Value, E> {
                Ok(UniqueJsonKeys)
            }

            fn visit_unit<E: de::Error>(self) -> Result<Self::Value, E> {
                Ok(UniqueJsonKeys)
            }
        }
        deserializer.deserialize_any(UniqueVisitor)
    }
}

fn parse_number(value: &str, field: &str) -> Result<u64, EngineError> {
    if value.is_empty()
        || (value.len() > 1 && value.starts_with('0'))
        || !value.bytes().all(|b| b.is_ascii_digit())
    {
        return Err(EngineError::invalid_input(format!(
            "FASTA {field} must be a canonical unsigned decimal integer"
        )));
    }
    value
        .parse()
        .map_err(|_| EngineError::invalid_input(format!("FASTA {field} integer overflow")))
}

fn parse_header(header: &str, reference: &TssReference) -> Result<TssRecord, EngineError> {
    if header.len() > MAX_HEADER_BYTES {
        return Err(EngineError::invalid_input("FASTA header exceeds 64 KiB"));
    }
    let mut fields = BTreeMap::new();
    for pair in header.split('|') {
        let (key, value) = pair
            .split_once('=')
            .ok_or_else(|| EngineError::invalid_input("FASTA headers require key=value fields"))?;
        if !HEADER_KEYS.contains(&key) || value.is_empty() || value.contains('=') {
            return Err(EngineError::invalid_input(format!(
                "Unknown, empty or malformed FASTA header field {key}"
            )));
        }
        if fields.insert(key, value).is_some() {
            return Err(EngineError::invalid_input(format!(
                "Duplicate FASTA header key {key}"
            )));
        }
    }
    if fields.len() != HEADER_KEYS.len() {
        return Err(EngineError::invalid_input(
            "FASTA header is missing required bundle fields",
        ));
    }
    if fields["orientation"] != ORIENTATION || fields["assembly"] != reference.assembly {
        return Err(EngineError::invalid_input(
            "FASTA assembly/orientation must match the manifest and transcript_5prime_to_3prime convention",
        ));
    }
    let (start, end) = fields["genomic_1based"].split_once("..").ok_or_else(|| {
        EngineError::invalid_input("FASTA genomic_1based must have START..END geometry")
    })?;
    let (upstream, downstream) = fields["window"]
        .split_once("..")
        .and_then(|(up, down)| Some((up.strip_prefix('-')?, down.strip_prefix('+')?)))
        .ok_or_else(|| EngineError::invalid_input("FASTA window must have -U..+D geometry"))?;
    let strand = match fields["strand"] {
        "+" => TssStrand::Plus,
        "-" => TssStrand::Minus,
        _ => return Err(EngineError::invalid_input("FASTA strand must be + or -")),
    };
    let record = TssRecord {
        promoter_id: fields["promoter_id"].into(),
        gene_id: fields["gene_id"].into(),
        gene_symbol: fields["gene_symbol"].into(),
        geometry: TssGeometry {
            chromosome: fields["chromosome"].into(),
            strand,
            tss_1based: parse_number(fields["tss_1based"], "tss_1based")?,
            start_1based: parse_number(start, "genomic start")?,
            end_1based: parse_number(end, "genomic end")?,
            upstream_bp: usize::try_from(parse_number(upstream, "upstream extent")?)
                .map_err(|_| EngineError::invalid_input("FASTA upstream extent overflow"))?,
            downstream_bp: usize::try_from(parse_number(downstream, "downstream extent")?)
                .map_err(|_| EngineError::invalid_input("FASTA downstream extent overflow"))?,
        },
        transcripts: fields["transcripts"].split(',').map(String::from).collect(),
        sequence_sha256: fields["sequence_sha256"].into(),
    };
    validate_record(&record)?;
    Ok(record)
}

fn parse_fasta(
    bytes: &[u8],
    name: &str,
    reference: &TssReference,
    expected: &BTreeMap<&str, &TssRecord>,
    sequences: &mut BTreeMap<String, String>,
) -> Result<(), EngineError> {
    parse_fasta_with_header(bytes, name, reference, expected, sequences, parse_header)
}

fn parse_fasta_with_header(
    bytes: &[u8],
    name: &str,
    reference: &TssReference,
    expected: &BTreeMap<&str, &TssRecord>,
    sequences: &mut BTreeMap<String, String>,
    header_parser: fn(&str, &TssReference) -> Result<TssRecord, EngineError>,
) -> Result<(), EngineError> {
    let text = std::str::from_utf8(bytes)
        .map_err(|_| EngineError::invalid_input(format!("FASTA {name} must be UTF-8 text")))?;
    let mut current: Option<(TssRecord, String)> = None;
    let mut record_count = 0;
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some((record, sequence)) = current.take() {
                finish_record(record, sequence, sequences)?;
            }
            let record = header_parser(header, reference)?;
            let manifest_record = expected.get(record.promoter_id.as_str()).ok_or_else(|| {
                EngineError::invalid_input(format!(
                    "FASTA promoter {} is absent from the manifest",
                    record.promoter_id
                ))
            })?;
            let transcripts: BTreeSet<_> = record.transcripts.iter().collect();
            if record.gene_id != manifest_record.gene_id
                || record.gene_symbol != manifest_record.gene_symbol
                || record.geometry != manifest_record.geometry
                || record.sequence_sha256 != manifest_record.sequence_sha256
                || transcripts != manifest_record.transcripts.iter().collect()
            {
                return Err(EngineError::invalid_input(format!(
                    "FASTA/manifest metadata mismatch for promoter {}",
                    record.promoter_id
                )));
            }
            record_count += 1;
            if record_count > expected.len() {
                return Err(EngineError::invalid_input(
                    "FASTA record count exceeds the manifest",
                ));
            }
            current = Some((record, String::new()));
        } else {
            for byte in line.bytes().filter(|byte| !byte.is_ascii_whitespace()) {
                let Some((record, sequence)) = current.as_mut() else {
                    return Err(EngineError::invalid_input(
                        "FASTA sequence content appears before its first header",
                    ));
                };
                let base = byte.to_ascii_uppercase();
                if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
                    return Err(EngineError::invalid_input(format!(
                        "FASTA sequence for {} contains a base outside A/C/G/T/N",
                        record.promoter_id
                    )));
                }
                if sequence.len() >= record.geometry.length().unwrap_or(0) {
                    return Err(EngineError::invalid_input(format!(
                        "FASTA sequence length exceeds geometry for {}",
                        record.promoter_id
                    )));
                }
                sequence.push(base as char);
            }
        }
    }
    if let Some((record, sequence)) = current {
        finish_record(record, sequence, sequences)?;
    }
    if record_count == 0 {
        return Err(EngineError::invalid_input(format!(
            "FASTA member {name} has no records"
        )));
    }
    Ok(())
}

fn finish_record(
    record: TssRecord,
    sequence: String,
    sequences: &mut BTreeMap<String, String>,
) -> Result<(), EngineError> {
    if Some(sequence.len()) != record.geometry.length() {
        return Err(EngineError::invalid_input(format!(
            "FASTA sequence length does not match geometry for {}",
            record.promoter_id
        )));
    }
    if sha256_hex_bytes(sequence.as_bytes()) != record.sequence_sha256 {
        return Err(EngineError::invalid_input(format!(
            "Normalized sequence SHA-256 mismatch for {}",
            record.promoter_id
        )));
    }
    if sequences
        .insert(record.promoter_id.clone(), sequence)
        .is_some()
    {
        return Err(EngineError::invalid_input(format!(
            "Duplicate FASTA promoter ID {}",
            record.promoter_id
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use gentle_protocol::ErrorCode;
    use gentle_protocol::tss_profiles::TssSelectedRecord;
    use tempfile::TempDir;

    const FIXTURE_MEMBERS: [&str; 4] = ["manifest.json", "minus.fa", "plus.fa", "selection.json"];

    fn fixture_root() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("test_files/fixtures/tss_profiles")
    }

    fn request(root: &Path) -> ComputeTssProfilesRequest {
        ComputeTssProfilesRequest {
            manifest: root.join("manifest.json").to_string_lossy().into_owned(),
            panel: "unused-by-the-bundle-reader.json".into(),
            fasta: vec![],
            selection: None,
            expected_genome_id: "synthetic-genome-v1".into(),
            expected_assembly: Some("synthetic-assembly-v1".into()),
            expected_annotation_release: Some("synthetic-annotation-v1".into()),
            expected_dataset_id: None,
        }
    }

    fn copied_fixture() -> TempDir {
        let temp = tempfile::tempdir().unwrap();
        for name in FIXTURE_MEMBERS.into_iter().chain([CHECKSUMS_NAME]) {
            fs::copy(fixture_root().join(name), temp.path().join(name)).unwrap();
        }
        temp
    }

    fn load_manifest(root: &Path) -> TssBundleManifest {
        serde_json::from_slice(&fs::read(root.join("manifest.json")).unwrap()).unwrap()
    }

    fn save_manifest(root: &Path, manifest: &TssBundleManifest) {
        let text = serde_json::to_string_pretty(manifest).unwrap() + "\n";
        fs::write(root.join("manifest.json"), text).unwrap();
    }

    fn refresh_checksums(root: &Path) {
        let sums: String = FIXTURE_MEMBERS
            .iter()
            .map(|name| {
                format!(
                    "{}  {name}\n",
                    sha256_hex_bytes(&fs::read(root.join(name)).unwrap())
                )
            })
            .collect();
        fs::write(root.join(CHECKSUMS_NAME), sums).unwrap();
    }

    fn rebind(root: &Path) {
        let mut manifest = load_manifest(root);
        for (name, digest) in &mut manifest.fasta_files {
            *digest = sha256_hex_bytes(&fs::read(root.join(name)).unwrap());
        }
        save_manifest(root, &manifest);
        refresh_checksums(root);
    }

    fn replace_plus(root: &Path, before: &str, after: &str) {
        let path = root.join("plus.fa");
        let original = fs::read_to_string(&path).unwrap();
        assert!(original.contains(before), "missing replacement {before}");
        fs::write(path, original.replace(before, after)).unwrap();
    }

    fn invalid(request: &ComputeTssProfilesRequest, diagnostic: &str) {
        let error = read_bundle(request).unwrap_err();
        assert_eq!(error.code, ErrorCode::InvalidInput);
        assert!(
            error.message.contains(diagnostic),
            "expected {diagnostic:?}, received {error:?}"
        );
    }

    #[test]
    fn synthetic_bundle_preserves_minus_orientation_memberships_and_explicit_selection() {
        let root = fixture_root();
        let mut input = request(&root);
        let unselected = read_bundle(&input).unwrap();
        assert!(unselected.records.iter().all(|(_, _, selected)| !selected));
        input.selection = Some("selection.json".into());
        let bundle = read_bundle(&input).unwrap();
        assert_eq!(bundle.reference.genome_id, "synthetic-genome-v1");
        assert_eq!(bundle.records.len(), 2);
        let (plus, plus_sequence, plus_selected) = &bundle.records[0];
        assert_eq!(plus.promoter_id, "synthetic-plus");
        assert_eq!(plus_sequence, "ACGTNACGTAC");
        assert!(!plus_selected);
        assert_eq!(
            plus.transcripts,
            ["synthetic.plus.tx2", "synthetic.plus.tx1"]
        );
        assert_eq!(plus.geometry.genomic_at(0), Some(97));
        assert_eq!(plus.geometry.genomic_at(3), Some(100));
        assert_eq!(plus.geometry.genomic_at(10), Some(107));
        let (minus, minus_sequence, minus_selected) = &bundle.records[1];
        assert_eq!(minus.promoter_id, "synthetic-minus");
        assert_eq!(minus_sequence, "TTACGACGTAA");
        assert!(*minus_selected);
        assert_eq!(minus.geometry.strand, TssStrand::Minus);
        assert_eq!(minus.geometry.genomic_at(0), Some(303));
        assert_eq!(minus.geometry.genomic_at(3), Some(300));
        assert_eq!(minus.geometry.genomic_at(10), Some(293));
        assert_eq!(minus.geometry.relative_at(0), Some(-3));
        assert_eq!(minus.geometry.relative_at(10), Some(7));
        assert_eq!(bundle.inputs.len(), 5);
        for input in &bundle.inputs {
            assert!(!Path::new(&input.name).is_absolute());
            assert_eq!(
                input.sha256,
                sha256_hex_bytes(&fs::read(root.join(&input.name)).unwrap())
            );
        }
        assert!(bundle.warnings[0].contains("not independently verified"));
        assert_eq!(bundle.source.schema, BUNDLE_SCHEMA);
        assert_eq!(
            bundle.source.manifest_sha256,
            sha256_hex_bytes(&fs::read(root.join("manifest.json")).unwrap())
        );
        assert!(bundle.source.source_revision.is_none());
        assert!(bundle.source.dataset_id.is_none());
        assert!(bundle.source.producer_sha256.is_none());
        assert!(bundle.selection_evidence.is_empty());
        let replay = read_bundle(&input).unwrap();
        assert_eq!(bundle.records, replay.records);
        assert_eq!(bundle.warnings, replay.warnings);
    }

    #[test]
    fn fasta_list_requires_exact_complete_paths_without_duplicates() {
        let root = fixture_root();
        let mut input = request(&root);
        input.fasta = vec!["plus.fa".into(), "minus.fa".into()];
        read_bundle(&input).unwrap();
        input.fasta = vec![
            root.join("minus.fa").to_string_lossy().into_owned(),
            root.join("plus.fa").to_string_lossy().into_owned(),
        ];
        read_bundle(&input).unwrap();
        for files in [
            vec!["plus.fa"],
            vec!["plus.fa", "plus.fa"],
            vec!["plus.fa", "minus.fa", "other.fa"],
            vec!["elsewhere/plus.fa", "minus.fa"],
            vec!["./plus.fa", "minus.fa"],
        ] {
            input.fasta = files.into_iter().map(String::from).collect();
            assert!(read_bundle(&input).is_err());
        }
    }

    #[test]
    fn byte_tampering_fails_for_manifest_fasta_and_selection() {
        for member in ["manifest.json", "plus.fa", "selection.json"] {
            let temp = copied_fixture();
            let path = temp.path().join(member);
            let mut bytes = fs::read(&path).unwrap();
            bytes.push(b'\n');
            fs::write(path, bytes).unwrap();
            let mut input = request(temp.path());
            input.selection = Some("selection.json".into());
            invalid(&input, "SHA256SUMS");
        }
    }

    #[test]
    fn sequence_digest_is_checked_independently_of_rebound_file_hashes() {
        let temp = copied_fixture();
        replace_plus(temp.path(), "acgtn ac", "tcgtn ac");
        rebind(temp.path());
        invalid(
            &request(temp.path()),
            "Normalized sequence SHA-256 mismatch",
        );
    }

    #[test]
    fn fastas_reject_invalid_alphabets_missing_bases_and_unheaded_data() {
        for sequence in [
            "ACGTUACGTAC",
            "ACGTRACGTAC",
            "ACGT-ACGTAC",
            "ACGTNACGTA",
            "ACGTNACGTACC",
            "ACGTN\u{a0}ACGTAC",
            "",
        ] {
            let temp = copied_fixture();
            replace_plus(temp.path(), "acgtn ac\ngtac", sequence);
            rebind(temp.path());
            assert!(read_bundle(&request(temp.path())).is_err(), "{sequence:?}");
        }
        let temp = copied_fixture();
        let original = fs::read_to_string(temp.path().join("plus.fa")).unwrap();
        fs::write(temp.path().join("plus.fa"), format!("ACGT\n{original}")).unwrap();
        rebind(temp.path());
        invalid(&request(temp.path()), "before its first header");
    }

    #[test]
    fn headers_reject_unknown_duplicate_missing_and_mismatched_fields() {
        for (before, after) in [
            ("gene_symbol=SYNPLUS", "gene_symbol=SYNPLUS|extra=x"),
            (
                "gene_symbol=SYNPLUS",
                "gene_symbol=SYNPLUS|gene_symbol=SYNPLUS",
            ),
            ("gene_symbol=SYNPLUS|", ""),
            ("gene_symbol=SYNPLUS", "gene_symbol=other"),
            ("gene_id=synthetic-gene-plus", "gene_id=other"),
            ("promoter_id=synthetic-plus", "promoter_id=unknown"),
            ("assembly=synthetic-assembly-v1", "assembly=other"),
            ("chromosome=synthetic_chr", "chromosome=other"),
            ("strand=+", "strand=-"),
            ("tss_1based=100", "tss_1based=0100"),
            ("tss_1based=100", "tss_1based=18446744073709551616"),
            ("genomic_1based=97..107", "genomic_1based=98..107"),
            ("window=-3..+7", "window=-4..+7"),
            ("window=-3..+7", "window=3..7"),
            (
                "orientation=transcript_5prime_to_3prime",
                "orientation=genomic",
            ),
            (
                "synthetic.plus.tx1,synthetic.plus.tx2",
                "synthetic.plus.tx1,synthetic.plus.tx1",
            ),
            (
                "synthetic.plus.tx1,synthetic.plus.tx2",
                "synthetic.plus.tx1,unknown",
            ),
            (
                "synthetic.plus.tx1,synthetic.plus.tx2",
                "synthetic.plus.tx1",
            ),
            ("sequence_sha256=4d9e", "sequence_sha256=4D9E"),
            ("sequence_sha256=4d9e", "sequence_sha256=sha256:4d9e"),
        ] {
            let temp = copied_fixture();
            replace_plus(temp.path(), before, after);
            rebind(temp.path());
            assert!(read_bundle(&request(temp.path())).is_err(), "{after}");
        }
    }

    #[test]
    fn duplicate_manifest_ids_and_physical_tsss_fail_without_deduplicating_memberships() {
        for duplicate_id in [true, false] {
            let temp = copied_fixture();
            let mut manifest = load_manifest(temp.path());
            let mut repeated = manifest.records[0].clone();
            if !duplicate_id {
                repeated.promoter_id = "another-promoter-at-same-tss".into();
                repeated.gene_id = "another-gene-at-same-tss".into();
            }
            manifest.records.push(repeated);
            save_manifest(temp.path(), &manifest);
            refresh_checksums(temp.path());
            invalid(
                &request(temp.path()),
                if duplicate_id {
                    "Duplicate promoter"
                } else {
                    "Duplicate physical TSS"
                },
            );
        }
        let temp = copied_fixture();
        let mut manifest = load_manifest(temp.path());
        manifest.records[0]
            .transcripts
            .push("synthetic.plus.tx1".into());
        save_manifest(temp.path(), &manifest);
        refresh_checksums(temp.path());
        invalid(&request(temp.path()), "Duplicate transcript");
    }

    #[test]
    fn duplicate_fasta_ids_and_missing_manifest_members_fail() {
        let temp = copied_fixture();
        let path = temp.path().join("plus.fa");
        let fasta = fs::read_to_string(&path).unwrap();
        fs::write(path, format!("{fasta}{fasta}")).unwrap();
        rebind(temp.path());
        invalid(&request(temp.path()), "Duplicate FASTA promoter");

        let temp = copied_fixture();
        let mut manifest = load_manifest(temp.path());
        manifest.fasta_files.remove("minus.fa");
        save_manifest(temp.path(), &manifest);
        refresh_checksums(temp.path());
        invalid(&request(temp.path()), "membership mismatch");
    }

    #[test]
    fn repeated_sequence_at_a_distinct_locus_is_retained_with_a_diagnostic() {
        let temp = copied_fixture();
        let mut manifest = load_manifest(temp.path());
        let mut repeated = manifest.records[0].clone();
        repeated.promoter_id = "synthetic-repeated-sequence".into();
        repeated.geometry.tss_1based += 400;
        repeated.geometry.start_1based += 400;
        repeated.geometry.end_1based += 400;
        repeated.transcripts = vec!["synthetic.repeat.tx1".into()];
        manifest.records.push(repeated);
        save_manifest(temp.path(), &manifest);
        let original = fs::read_to_string(temp.path().join("plus.fa")).unwrap();
        let copy = original
            .replace(
                "promoter_id=synthetic-plus",
                "promoter_id=synthetic-repeated-sequence",
            )
            .replace("tss_1based=100", "tss_1based=500")
            .replace("genomic_1based=97..107", "genomic_1based=497..507")
            .replace(
                "synthetic.plus.tx1,synthetic.plus.tx2",
                "synthetic.repeat.tx1",
            );
        fs::write(temp.path().join("plus.fa"), format!("{original}{copy}")).unwrap();
        rebind(temp.path());
        let bundle = read_bundle(&request(temp.path())).unwrap();
        assert_eq!(bundle.records.len(), 3);
        assert_eq!(bundle.records[0].1, bundle.records[2].1);
        assert_ne!(bundle.records[0].0.geometry, bundle.records[2].0.geometry);
        assert!(
            bundle
                .warnings
                .iter()
                .any(|warning| warning.contains("Repeated sequence"))
        );
    }

    #[test]
    fn reference_expectations_are_independent_exact_checks() {
        let root = fixture_root();
        for field in 0..3 {
            let mut input = request(&root);
            match field {
                0 => input.expected_genome_id = "synthetic-assembly-v1".into(),
                1 => input.expected_assembly = Some("other".into()),
                _ => input.expected_annotation_release = Some("other".into()),
            }
            invalid(&input, "reference does not exactly match");
        }
        let mut input = request(&root);
        input.expected_assembly = None;
        input.expected_annotation_release = None;
        read_bundle(&input).unwrap();
        input.expected_dataset_id = Some("not-declared-by-canonical-schema".into());
        invalid(&input, "does not declare a dataset_id");
    }

    #[test]
    fn selection_requires_known_noncontradictory_unique_references() {
        for variant in 0..5 {
            let temp = copied_fixture();
            let path = temp.path().join("selection.json");
            let mut selection: TssSelection =
                serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
            match variant {
                0 => selection.selected[0].promoter_id = "missing".into(),
                1 => selection.selected[0].gene_id = "synthetic-gene-plus".into(),
                2 => selection.reference.annotation_release = Some("other".into()),
                3 => selection.selected.push(TssSelectedRecord {
                    promoter_id: selection.selected[0].promoter_id.clone(),
                    gene_id: selection.selected[0].gene_id.clone(),
                }),
                _ => selection.schema = "legacy.unspecified".into(),
            }
            fs::write(path, serde_json::to_vec(&selection).unwrap()).unwrap();
            refresh_checksums(temp.path());
            let mut input = request(temp.path());
            input.selection = Some("selection.json".into());
            assert!(read_bundle(&input).is_err());
        }
        let temp = copied_fixture();
        let path = temp.path().join(CHECKSUMS_NAME);
        let sums = fs::read_to_string(&path).unwrap();
        let sums: String = sums
            .lines()
            .filter(|line| !line.ends_with("selection.json"))
            .map(|line| format!("{line}\n"))
            .collect();
        fs::write(path, sums).unwrap();
        let mut input = request(temp.path());
        input.selection = Some("selection.json".into());
        invalid(&input, "listed in SHA256SUMS");
    }

    #[test]
    fn checksum_inventory_rejects_missing_duplicate_malformed_and_unsafe_members() {
        let valid = format!("{}  plus.fa\n", "a".repeat(64));
        for sums in [
            format!("{valid}{valid}"),
            format!("{}  plus.fa\n", "A".repeat(64)),
            format!("sha256:{}  plus.fa\n", "a".repeat(64)),
            format!("{} plus.fa\n", "a".repeat(64)),
            format!("{}  ../escape\n", "a".repeat(64)),
            format!("{}  /absolute\n", "a".repeat(64)),
            format!("{}  C:\\escape\n", "a".repeat(64)),
            format!("{}  ./plus.fa\n", "a".repeat(64)),
            format!("{}  SHA256SUMS\n", "a".repeat(64)),
        ] {
            assert!(parse_checksums(sums.as_bytes()).is_err(), "{sums}");
        }
        parse_checksums(format!("{} *plus.fa\r\n", "a".repeat(64)).as_bytes()).unwrap();
        let temp = copied_fixture();
        fs::write(temp.path().join(CHECKSUMS_NAME), valid).unwrap();
        invalid(
            &request(temp.path()),
            "missing or mismatched checksum for manifest.json",
        );
    }

    #[test]
    fn manifest_rejects_unknown_fields_duplicate_json_keys_and_unsafe_links() {
        let temp = copied_fixture();
        let path = temp.path().join("manifest.json");
        let original = fs::read_to_string(&path).unwrap();
        fs::write(&path, original.replacen('{', "{\"unknown\":true,", 1)).unwrap();
        refresh_checksums(temp.path());
        invalid(&request(temp.path()), "Invalid bundle manifest schema");
        let duplicate = original.replace(
            "\"fasta_files\": {",
            &format!("\"fasta_files\": {{\"plus.fa\":\"{}\",", "a".repeat(64)),
        );
        fs::write(path, duplicate).unwrap();
        refresh_checksums(temp.path());
        invalid(&request(temp.path()), "Invalid bundle manifest JSON");

        for name in [
            "../escape.fa",
            "./plus.fa",
            "nested/../plus.fa",
            "/absolute.fa",
            "a\\b.fa",
        ] {
            let temp = copied_fixture();
            let mut manifest = load_manifest(temp.path());
            let digest = manifest.fasta_files.remove("plus.fa").unwrap();
            manifest.fasta_files.insert(name.into(), digest);
            save_manifest(temp.path(), &manifest);
            refresh_checksums(temp.path());
            invalid(&request(temp.path()), "portable relative filename");
        }
    }

    #[test]
    fn byte_and_record_limits_fail_without_large_allocations() {
        let temp = copied_fixture();
        let large = File::create(temp.path().join("manifest.json")).unwrap();
        large.set_len(MAX_MANIFEST_BYTES + 1).unwrap();
        invalid(&request(temp.path()), "byte limit");
        let small = temp.path().join("bounded.txt");
        fs::write(&small, b"ACGT").unwrap();
        assert!(read_bounded(&small, 3, &mut 0).is_err());
        assert!(read_bounded(&small, 4, &mut (MAX_TOTAL_BYTES - 3)).is_err());
        let temp = copied_fixture();
        let mut manifest = load_manifest(temp.path());
        let record = &mut manifest.records[0];
        record.geometry.downstream_bp = MAX_SEQUENCE_BASES;
        record.geometry.end_1based = record.geometry.tss_1based + MAX_SEQUENCE_BASES as u64;
        save_manifest(temp.path(), &manifest);
        refresh_checksums(temp.path());
        invalid(&request(temp.path()), "base window limit");
        let temp = copied_fixture();
        let mut manifest = load_manifest(temp.path());
        manifest.records = vec![manifest.records[0].clone(); 10_001];
        let error = validate_manifest(&manifest).unwrap_err();
        assert!(error.message.contains("1..=10000"));
    }

    #[cfg(unix)]
    #[test]
    fn tss_regular_input_rejects_an_unopened_fifo_without_waiting_for_a_writer() {
        use std::os::unix::ffi::OsStrExt;
        let temp = tempfile::tempdir().unwrap();
        let fifo = temp.path().join("panel.json");
        let path = std::ffi::CString::new(fifo.as_os_str().as_bytes()).unwrap();
        // The NUL-terminated path stays alive during mkfifo; only this temporary
        // fixture directory is affected, and no writer is opened on the pipe.
        assert_eq!(unsafe { libc::mkfifo(path.as_ptr(), 0o600) }, 0);
        let (send, receive) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let _ = send.send(open_regular_input(&fifo).map(|_| ()).map_err(|e| e.message));
        });
        let error = receive
            .recv_timeout(std::time::Duration::from_secs(5))
            .expect("reject the pipe without blocking on an absent writer")
            .unwrap_err();
        assert!(error.contains("regular file"), "{error}");
        assert!(open_regular_input(temp.path()).is_err());
    }

    #[cfg(unix)]
    #[test]
    fn symlink_escapes_and_nonregular_members_are_rejected() {
        use std::os::unix::fs::symlink;

        for name in ["manifest.json", "plus.fa", "selection.json", CHECKSUMS_NAME] {
            let temp = copied_fixture();
            let outside = tempfile::tempdir().unwrap();
            fs::copy(temp.path().join(name), outside.path().join(name)).unwrap();
            fs::remove_file(temp.path().join(name)).unwrap();
            symlink(outside.path().join(name), temp.path().join(name)).unwrap();
            invalid(&request(temp.path()), "escapes its directory");
        }
        let temp = copied_fixture();
        fs::remove_file(temp.path().join("plus.fa")).unwrap();
        fs::create_dir(temp.path().join("plus.fa")).unwrap();
        invalid(&request(temp.path()), "regular file");
    }
}
