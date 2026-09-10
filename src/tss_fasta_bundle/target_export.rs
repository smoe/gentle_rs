//! Explicit adapter for the target-export v1 files inspected at source revision
//! `a106cbbd5223f8c4be55a7d846b8416bc4b3ae33`. No historical format is guessed.
//!
//! The manifest hashes SHA256SUMS, which hashes the original FASTAs. Its header
//! grammar uses a bare gene symbol, START-END, and minusU_plusD. Only in-memory
//! header metadata is adapted; source bytes and their digests remain untouched.

use super::*;

pub(super) const SCHEMA: &str = "gentle.target_tss_fasta_export.v1";
const INTEGRATED_SELECTION_SCHEMA: &str = "gentle.regulatory_region_comparison_sequences.v1";

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct TargetManifest {
    schema: String,
    assembly_id: String,
    upstream_bp: usize,
    downstream_bp: usize,
    sequence_orientation: String,
    files: Vec<TargetFile>,
    sha256sums_sha256: String,
    source: TargetSource,
    source_revision: String,
    producer_sha256: String,
    record_policy: String,
    non_claims: Vec<String>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct TargetSource {
    dataset_id: String,
    genome_id: String,
    promoter_fasta_sha256: String,
    promoter_transcripts_sha256: String,
    promoter_windows_sha256: String,
    promoterome_receipt_sha256: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct TargetFile {
    filename: String,
    gene_id: String,
    gene_symbol: String,
    record_count: usize,
    records: Vec<TargetRecord>,
    sha256: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct TargetRecord {
    promoter_id: String,
    gene_id: String,
    chromosome: String,
    strand: TssStrand,
    tss_1based: u64,
    genomic_start_1based: u64,
    genomic_end_1based: u64,
    sequence_length_bp: usize,
    sequence_sha256: String,
    transcript_ids: Vec<String>,
}

pub(super) fn read_target_bundle(
    request: &ComputeTssProfilesRequest,
    mut files: BundleFiles,
) -> Result<VerifiedTssBundle, EngineError> {
    let raw: TargetManifest = parse_json(&files.manifest_bytes, "target-export manifest")?;
    if raw.schema != SCHEMA || raw.sequence_orientation != ORIENTATION {
        return Err(EngineError::invalid_input(
            "Target-export schema/orientation must match its explicit v1 contract",
        ));
    }
    if request.expected_assembly.is_none() || request.expected_dataset_id.is_none() {
        return Err(EngineError::invalid_input(
            "Target-export bundles require expected_assembly and expected_dataset_id",
        ));
    }
    text_field(&raw.source.dataset_id, "source.dataset_id")?;
    if request.expected_dataset_id.as_ref() != Some(&raw.source.dataset_id) {
        return Err(EngineError::invalid_input(
            "Target-export source.dataset_id does not exactly match expected_dataset_id",
        ));
    }
    let reference = TssReference {
        genome_id: raw.source.genome_id.clone(),
        assembly: raw.assembly_id.clone(),
        // This version declares no separate annotation release.
        annotation_release: None,
    };
    check_reference_expectations(request, &reference)?;
    validate_sha256(&raw.sha256sums_sha256)?;
    if sha256_hex_bytes(&files.sums_bytes) != raw.sha256sums_sha256 {
        return Err(EngineError::invalid_input(
            "Target-export manifest sha256sums_sha256 does not match SHA256SUMS bytes",
        ));
    }
    validate_sha256(&raw.producer_sha256)?;
    text_field(&raw.source_revision, "source_revision")?;
    text_field(&raw.record_policy, "record_policy")?;
    for claim in &raw.non_claims {
        text_field(claim, "non_claims")?;
    }

    let mut inputs = vec![
        binding(
            "bundle_manifest",
            &files.manifest_name,
            &files.manifest_bytes,
        ),
        binding("bundle_checksums", CHECKSUMS_NAME, &files.sums_bytes),
    ];
    // These are declarations in the hash-bound manifest, not files independently
    // read or verified by this adapter. Keep that distinction in the role.
    for (name, digest) in [
        ("promoter_fasta_sha256", &raw.source.promoter_fasta_sha256),
        (
            "promoter_transcripts_sha256",
            &raw.source.promoter_transcripts_sha256,
        ),
        (
            "promoter_windows_sha256",
            &raw.source.promoter_windows_sha256,
        ),
        (
            "promoterome_receipt_sha256",
            &raw.source.promoterome_receipt_sha256,
        ),
    ] {
        validate_sha256(digest)?;
        inputs.push(TssInputBinding {
            role: "declared_source".into(),
            name: name.into(),
            sha256: digest.clone(),
        });
    }

    let mut manifest = TssBundleManifest {
        schema: BUNDLE_SCHEMA.into(),
        reference,
        fasta_files: BTreeMap::new(),
        records: Vec::new(),
    };
    let mut members = BTreeMap::new();
    for file in &raw.files {
        validate_relative_name(&file.filename)?;
        if file.filename == files.manifest_name || file.filename == CHECKSUMS_NAME {
            return Err(EngineError::invalid_input(
                "Target-export FASTA cannot reuse manifest/checksum filenames",
            ));
        }
        if manifest
            .fasta_files
            .insert(file.filename.clone(), file.sha256.clone())
            .is_some()
        {
            return Err(EngineError::invalid_input(format!(
                "Duplicate target-export FASTA path {}",
                file.filename
            )));
        }
        if file.records.is_empty() || file.record_count != file.records.len() {
            return Err(EngineError::invalid_input(format!(
                "Target-export record_count mismatch or empty FASTA {}",
                file.filename
            )));
        }
        let first = manifest.records.len();
        for row in &file.records {
            if row.gene_id != file.gene_id {
                return Err(EngineError::invalid_input(format!(
                    "Target-export gene_id mismatch for promoter {} in {}",
                    row.promoter_id, file.filename
                )));
            }
            let record = TssRecord {
                promoter_id: row.promoter_id.clone(),
                gene_id: row.gene_id.clone(),
                gene_symbol: file.gene_symbol.clone(),
                geometry: TssGeometry {
                    chromosome: row.chromosome.clone(),
                    strand: row.strand,
                    tss_1based: row.tss_1based,
                    start_1based: row.genomic_start_1based,
                    end_1based: row.genomic_end_1based,
                    upstream_bp: raw.upstream_bp,
                    downstream_bp: raw.downstream_bp,
                },
                transcripts: row.transcript_ids.clone(),
                sequence_sha256: row.sequence_sha256.clone(),
            };
            if Some(row.sequence_length_bp) != record.geometry.length()
                || row.sequence_length_bp > MAX_SEQUENCE_BASES
            {
                return Err(EngineError::invalid_input(format!(
                    "Target-export sequence_length_bp does not match bounded geometry for promoter {}",
                    row.promoter_id
                )));
            }
            validate_record(&record).map_err(|error| {
                EngineError::invalid_input(format!(
                    "Target-export promoter {}: {}",
                    row.promoter_id, error.message
                ))
                .with_cause(error)
            })?;
            manifest.records.push(record);
        }
        members.insert(file.filename.as_str(), first..manifest.records.len());
    }
    validate_manifest(&manifest)?;
    check_requested_fastas(
        request,
        &files.root,
        &files.declared_parent,
        &manifest.fasta_files,
    )?;
    // In this schema the list contains exactly files[], never the manifest: a
    // manifest entry would make the declared checksum-list digest circular.
    if files.sums != manifest.fasta_files {
        return Err(EngineError::invalid_input(
            "Target-export SHA256SUMS must exactly match manifest FASTA paths and digests (no circular manifest entry)",
        ));
    }

    let mut resolved_paths = BTreeSet::from([files.manifest_path, files.sums_path]);
    let mut sequences = BTreeMap::new();
    for (name, range) in members {
        let path = contained_file(&files.root, name)?;
        if !resolved_paths.insert(path.clone()) {
            return Err(EngineError::invalid_input(format!(
                "Target-export bundle paths alias the same file: {name}"
            )));
        }
        let bytes = read_bounded(&path, MAX_MEMBER_BYTES, &mut files.total_bytes)?;
        verify_checksum(&files.sums, name, &bytes)?;
        let expected: BTreeMap<_, _> = manifest.records[range]
            .iter()
            .map(|row| (row.promoter_id.as_str(), row))
            .collect();
        let before = sequences.len();
        parse_fasta_with_header(
            &bytes,
            name,
            &manifest.reference,
            &expected,
            &mut sequences,
            parse_target_header,
        )?;
        if sequences.len() - before != expected.len() {
            return Err(EngineError::invalid_input(format!(
                "Target-export FASTA/manifest record membership mismatch in {name}"
            )));
        }
        inputs.push(binding("fasta", name, &bytes));
    }

    let selection_evidence = if let Some(value) = &request.selection {
        let (path, name) = selection_file(&files.root, &files.declared_parent, value)?;
        if !resolved_paths.insert(path.clone()) {
            return Err(EngineError::invalid_input(
                "Target-export selection must be distinct from manifest, checksums and FASTA files",
            ));
        }
        let bytes = read_bounded(&path, MAX_MANIFEST_BYTES, &mut files.total_bytes)?;
        let evidence = read_selection(&bytes, &manifest)?;
        inputs.push(binding("selection", &name, &bytes));
        evidence
    } else {
        BTreeMap::new()
    };
    let selected = selection_evidence.keys().cloned().collect();
    let mut warnings = raw.non_claims;
    warnings.push("Upstream declared_source digests, source revision and producer identity are preserved from the manifest; their source artifacts were not independently verified.".into());
    finish_bundle(
        manifest.records,
        sequences,
        &selected,
        VerifiedTssBundle {
            reference: manifest.reference,
            source: TssBundleSource {
                schema: raw.schema,
                manifest_sha256: sha256_hex_bytes(&files.manifest_bytes),
                source_revision: Some(raw.source_revision),
                dataset_id: Some(raw.source.dataset_id),
                producer_sha256: Some(raw.producer_sha256),
            },
            records: Vec::new(),
            inputs,
            warnings,
            selection_evidence,
        },
    )
}

fn parse_target_header(header: &str, reference: &TssReference) -> Result<TssRecord, EngineError> {
    if header.len() > MAX_HEADER_BYTES {
        return Err(EngineError::invalid_input("FASTA header exceeds 64 KiB"));
    }
    let (symbol, rest) = header.split_once('|').ok_or_else(|| {
        EngineError::invalid_input(
            "Target-export FASTA header requires a bare gene symbol followed by fields",
        )
    })?;
    let mut canonical = format!("gene_symbol={symbol}");
    for field in rest.split('|') {
        let (key, value) = field.split_once('=').ok_or_else(|| {
            EngineError::invalid_input("Target-export FASTA headers require key=value fields")
        })?;
        canonical.push('|');
        canonical.push_str(key);
        canonical.push('=');
        match key {
            "genomic_1based" => {
                let (start, end) = value.split_once('-').ok_or_else(|| {
                    EngineError::invalid_input(
                        "Target-export genomic_1based must have START-END geometry",
                    )
                })?;
                canonical.push_str(&format!("{start}..{end}"));
            }
            "window" => {
                let (up, down) = value
                    .strip_prefix("minus")
                    .and_then(|v| v.split_once("_plus"))
                    .ok_or_else(|| {
                        EngineError::invalid_input(
                            "Target-export window must have minusU_plusD geometry",
                        )
                    })?;
                canonical.push_str(&format!("-{up}..+{down}"));
            }
            _ => canonical.push_str(value),
        }
    }
    parse_header(&canonical, reference)
}

fn selection_file(
    root: &Path,
    declared_parent: &Path,
    value: &str,
) -> Result<(PathBuf, String), EngineError> {
    let path = Path::new(value);
    if !path.is_absolute() || path.starts_with(root) || path.starts_with(declared_parent) {
        let name = input_name(root, declared_parent, value)?;
        return Ok((contained_file(root, &name)?, name));
    }
    // Outside-bundle access is authorized only by this explicit absolute input,
    // never a filename embedded in the manifest or selection document.
    // Components normalize interior '.', so reject its original spelling too.
    if value
        .split(std::path::is_separator)
        .any(|part| matches!(part, "." | ".."))
    {
        return Err(EngineError::invalid_input(
            "Selection path must not contain traversal components",
        ));
    }
    for component in path.components() {
        match component {
            std::path::Component::Prefix(_) | std::path::Component::RootDir => {}
            std::path::Component::Normal(name) => validate_relative_name(
                name.to_str()
                    .ok_or_else(|| EngineError::invalid_input("Selection path must be UTF-8"))?,
            )?,
            std::path::Component::ParentDir | std::path::Component::CurDir => {
                return Err(EngineError::invalid_input(
                    "Selection path must not contain traversal components",
                ));
            }
        }
    }
    let name = path
        .file_name()
        .and_then(|v| v.to_str())
        .ok_or_else(|| EngineError::invalid_input("Selection needs a UTF-8 filename"))?
        .to_owned();
    validate_relative_name(&name)?;
    let path = fs::canonicalize(path).map_err(|error| {
        EngineError::invalid_input("Cannot resolve explicitly supplied selection file")
            .with_cause(error)
    })?;
    if !fs::metadata(&path)
        .map_err(|error| {
            EngineError::invalid_input("Cannot inspect selection file").with_cause(error)
        })?
        .is_file()
    {
        return Err(EngineError::invalid_input(
            "Selection must be a regular file",
        ));
    }
    Ok((path, name))
}

// Only identity and descriptive selection fields are interpreted. Other fields,
// notably evidence-window extents, length and digest, describe a different
// sequence and must never participate in the display-window join. Their original
// representation, including numeric evidence, is bound by the whole-file hash.
#[derive(Deserialize)]
struct IntegratedSelection {
    schema: String,
    regions: Vec<SelectedRegion>,
}

#[derive(Deserialize)]
struct SelectedRegion {
    promoterome_id: String,
    input_kind: String,
    transcript_ids: Vec<String>,
    genome_extraction: SelectedExtraction,
    cutrun_support: CutrunSupport,
}

#[derive(Deserialize)]
struct SelectedExtraction {
    gene_id: String,
    gene_name: String,
    genome_id: String,
    chromosome: String,
    strand: TssStrand,
    tss_1based: u64,
    transcript_ids: Vec<String>,
}

#[derive(Deserialize)]
struct CutrunSupport {
    criterion: String,
    factor: Option<String>,
    non_claim: String,
}

fn read_selection(
    bytes: &[u8],
    manifest: &TssBundleManifest,
) -> Result<BTreeMap<String, TssSelectionEvidence>, EngineError> {
    let selection: IntegratedSelection = parse_json(bytes, "integrated TSS selection")?;
    if selection.schema != INTEGRATED_SELECTION_SCHEMA {
        return Err(EngineError::invalid_input(format!(
            "Expected {INTEGRATED_SELECTION_SCHEMA} selection for target-export bundle"
        )));
    }
    let expected: BTreeMap<_, _> = manifest
        .records
        .iter()
        .map(|row| (row.promoter_id.as_str(), row))
        .collect();
    if selection.regions.len() > expected.len() {
        return Err(EngineError::invalid_input(
            "Selected promoter count exceeds the bundle",
        ));
    }
    let mut result = BTreeMap::new();
    for region in selection.regions {
        let id = &region.promoterome_id;
        let record = expected.get(id.as_str()).ok_or_else(|| {
            EngineError::invalid_input(format!(
                "Selected promoterome_id {id} is absent from the bundle"
            ))
        })?;
        if result.contains_key(id) {
            return Err(EngineError::invalid_input(format!(
                "Duplicate selected promoterome_id {id}"
            )));
        }
        let extraction = region.genome_extraction;
        if region.input_kind != "transcript_promoter"
            || extraction.gene_id != record.gene_id
            || extraction.gene_name != record.gene_symbol
            || extraction.genome_id != manifest.reference.genome_id
            || extraction.chromosome != record.geometry.chromosome
            || extraction.strand != record.geometry.strand
            || extraction.tss_1based != record.geometry.tss_1based
        {
            return Err(EngineError::invalid_input(format!(
                "Selection gene/genome/TSS/strand identity mismatch for promoter {id}"
            )));
        }
        transcript_subset(
            id,
            &region.transcript_ids,
            &record.transcripts,
            "transcript_ids",
        )?;
        transcript_subset(
            id,
            &extraction.transcript_ids,
            &record.transcripts,
            "genome_extraction.transcript_ids",
        )?;
        let support = region.cutrun_support;
        text_field(&support.criterion, "cutrun_support.criterion")?;
        text_field(&support.non_claim, "cutrun_support.non_claim")?;
        if let Some(factor) = &support.factor {
            text_field(factor, "cutrun_support.factor")?;
        }
        let factor_label = support
            .factor
            .as_ref()
            .map(|s| format!("{s} "))
            .unwrap_or_default();
        result.insert(region.promoterome_id, TssSelectionEvidence {
            label: format!("Selected in integrated report \u{2014} {factor_label}CUT&RUN-supported TSS window"),
            legend: format!("This TSS belonged to the integrated reporter panel and met the recorded descriptive CUT&RUN window criterion. This does not establish TSS usage, direct binding, or promoter activity. Source non-claim: {}", support.non_claim),
            criterion: support.criterion,
            factor: support.factor,
        });
    }
    Ok(result)
}

fn transcript_subset(
    id: &str,
    transcripts: &[String],
    manifest: &[String],
    field: &str,
) -> Result<(), EngineError> {
    if transcripts.is_empty() || transcripts.len() > 4_096 {
        return Err(EngineError::invalid_input(format!(
            "Selection {field} requires 1..=4096 memberships for promoter {id}"
        )));
    }
    let expected: BTreeSet<_> = manifest.iter().collect();
    let mut seen = BTreeSet::new();
    for transcript in transcripts {
        if !seen.insert(transcript) {
            return Err(EngineError::invalid_input(format!(
                "Duplicate selection transcript {transcript} in {field} for promoter {id}"
            )));
        }
        if !expected.contains(transcript) {
            return Err(EngineError::invalid_input(format!(
                "Selection {field} is not a transcript subset for promoter {id}: {transcript}"
            )));
        }
    }
    Ok(())
}

fn text_field(value: &str, field: &str) -> Result<(), EngineError> {
    if value.is_empty()
        || value.trim() != value
        || value.len() > 4_096
        || value.chars().any(char::is_control)
    {
        return Err(EngineError::invalid_input(format!(
            "{field} must be nonempty trimmed text without controls (at most 4096 bytes)"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
