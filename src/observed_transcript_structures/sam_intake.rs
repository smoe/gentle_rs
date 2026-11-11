//! Private verified SAM/dictionary/manifest intake with exact declared-unit accounting.
//!
//! No output publication, command dispatch, original-read identity or reference-byte
//! verification. Producer search, conversion and mask provenance remain declarations.

use super::*;
use crate::engine::{
    PrimerPayloadFile, REFERENCE_CONTIG_DICTIONARY_SCHEMA, ReferenceContigDictionary,
};
use std::fs;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

#[cfg(test)]
#[path = "sam_intake_tests.rs"]
pub(super) mod tests;

const MANIFEST_SCHEMA: &str = "gentle.observed_sam_query_units.v1";

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Conversion {
    DirectSam {
        producer: String,
        version: String,
        arguments: Vec<String>,
        selection_scope: String,
    },
    BamConversion {
        bam_sha256: String,
        producer: String,
        version: String,
        arguments: Vec<String>,
        selection_scope: String,
    },
    /// GENtle-executed indexed selection, distinct from a producer's declaration.
    ExecutedIndexedSelection {
        receipt: Box<super::indexed_locus::Receipt>,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaskWitness {
    pub coordinate_convention: String,
    pub unclipped_axis_length: usize,
    pub intervals: Vec<Interval>,
    pub source_sha256: String,
    pub derivation: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SearchDeclaration {
    pub search_complete: bool,
    pub views_complete: bool,
    pub limitations: Vec<String>,
    pub producer: String,
    pub version: String,
    pub arguments: Vec<String>,
    pub selection_scope: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RowWitness {
    pub ordinal: u64,
    pub qname: String,
    pub read_group: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Unit {
    pub unit_id: String,
    pub unclipped_axis_length: Option<usize>,
    pub rows: Vec<RowWitness>,
    pub masks: Option<MaskWitness>,
    #[serde(default)]
    pub ends: [Vec<EndConstraint>; 2],
    pub search: Option<SearchDeclaration>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Manifest {
    pub schema: String,
    pub sam_sha256: String,
    pub units: Vec<Unit>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Request {
    pub sam: PrimerPayloadFile,
    pub manifest: PrimerPayloadFile,
    pub dictionary: PrimerPayloadFile,
    pub conversion: Conversion,
    pub reference_id: String,
    pub reference: ReferenceBinding,
    pub reference_gene: ReferenceGene,
    pub comparison_genes: Vec<ReferenceGene>,
    // The public service must verify annotation files/receipts before supplying these values.
    pub annotation_sources: Vec<AnnotationSource>,
    pub version_policy: sam_text::VersionPolicy,
    pub unpadded_reference: bool,
    pub required_flank_bases: usize,
    pub limits: sam_text::Limits,
    pub max_units: usize,
    pub kernel_budgets: Budgets,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MissingCompanion {
    pub unit_id: String,
    pub row_ordinal: u64,
    pub companion: sam_text::Companion,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Receipt {
    pub sam_sha256: String,
    pub manifest_sha256: String,
    pub manifest_semantic_sha256: String,
    pub dictionary_sha256: String,
    pub source_eof_verified: bool,
    pub coordinate_authority: String,
    pub identity_scope: String,
    pub version_policy: sam_text::VersionPolicy,
    pub unpadded_reference: bool,
    pub declared_version: Option<String>,
    pub counts: sam_text::Counts,
    pub query_units: usize,
    pub query_unit_manifest: Manifest,
    pub missing_sa_companions: Vec<MissingCompanion>,
    pub conversion: Conversion,
    pub retained_peak_bytes: usize,
    pub work: usize,
}

pub(super) struct Prepared {
    pub input: sam::SamInput,
    pub receipt: Receipt,
}

fn fail(reason: &str) -> EngineError {
    bad(&format!("observed_sam_intake_{reason}"))
}
fn sha(value: &str) -> bool {
    value.strip_prefix("sha256:").is_some_and(|s| {
        s.len() == 64
            && s.bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    })
}
fn text(value: &str) -> bool {
    !value.trim().is_empty() && value.len() <= 4096 && !value.chars().any(char::is_control)
}
fn identifier(value: &str) -> bool {
    text(value) && !value.chars().any(char::is_whitespace) && value != "*"
}

fn provenance(producer: &str, version: &str, arguments: &[String], scope: &str) -> bool {
    text(producer)
        && text(version)
        && text(scope)
        && arguments.len() <= 1024
        && arguments
            .iter()
            .all(|s| s.len() <= 4096 && !s.contains('\0'))
}

/// Validate the spelling before Path components can erase unsafe dot/empty components.
pub(crate) fn path(raw: &str) -> Result<&Path, EngineError> {
    if raw.trim().is_empty() || raw.len() > 32768 || raw.chars().any(char::is_control) {
        return Err(fail("raw_path"));
    }
    let path = Path::new(raw);
    let prefix = match path.components().next() {
        Some(std::path::Component::Prefix(v)) => v.as_os_str().len(),
        _ => 0,
    };
    let separator = |c| c == '/' || (cfg!(windows) && c == '\\');
    let mut body = &raw[prefix..];
    if body.starts_with(separator) {
        body = &body[1..];
    }
    if (prefix > 0 && !path.has_root())
        || body.is_empty()
        || body.split(separator).any(|s| {
            s.is_empty()
                || s == "."
                || s == ".."
                || (cfg!(windows)
                    && (s.contains(':')
                        || s.ends_with(['.', ' '])
                        || matches!(
                            s.split('.')
                                .next()
                                .unwrap_or("")
                                .to_ascii_uppercase()
                                .as_str(),
                            "CON"
                                | "PRN"
                                | "AUX"
                                | "NUL"
                                | "COM1"
                                | "COM2"
                                | "COM3"
                                | "COM4"
                                | "COM5"
                                | "COM6"
                                | "COM7"
                                | "COM8"
                                | "COM9"
                                | "LPT1"
                                | "LPT2"
                                | "LPT3"
                                | "LPT4"
                                | "LPT5"
                                | "LPT6"
                                | "LPT7"
                                | "LPT8"
                                | "LPT9"
                        )))
        })
    {
        return Err(fail("raw_path"));
    }
    Ok(path)
}

struct Verified {
    file: fs::File,
    owner: same_file::Handle,
    expected: PrimerPayloadFile,
    modified: Option<std::time::SystemTime>,
    bytes: u64,
    hash: ring::digest::Context,
}
impl Verified {
    fn open(bound: &PrimerPayloadFile) -> Result<Self, EngineError> {
        let path = path(&bound.path)?;
        let before = fs::symlink_metadata(path).map_err(|e| fail(&format!("file: {e}")))?;
        if !sha(&bound.sha256)
            || !before.is_file()
            || before.file_type().is_symlink()
            || before.len() != bound.bytes
        {
            return Err(fail("file_binding"));
        }
        let file = fs::File::open(path).map_err(|e| fail(&format!("file: {e}")))?;
        let owner =
            same_file::Handle::from_file(file.try_clone().map_err(|e| fail(&e.to_string()))?)
                .map_err(|e| fail(&e.to_string()))?;
        let opened = file.metadata().map_err(|e| fail(&e.to_string()))?;
        if same_file::Handle::from_path(path).map_err(|e| fail(&e.to_string()))? != owner
            || opened.len() != before.len()
            || opened.modified().ok() != before.modified().ok()
        {
            return Err(fail("file_changed"));
        }
        Ok(Self {
            file,
            owner,
            expected: bound.clone(),
            modified: before.modified().ok(),
            bytes: 0,
            hash: ring::digest::Context::new(&ring::digest::SHA256),
        })
    }
    fn finish(self) -> Result<(), EngineError> {
        let after = fs::symlink_metadata(&self.expected.path).map_err(|e| fail(&e.to_string()))?;
        let opened = self.file.metadata().map_err(|e| fail(&e.to_string()))?;
        let hash = format!(
            "sha256:{}",
            self.hash
                .finish()
                .as_ref()
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect::<String>()
        );
        if !after.is_file()
            || after.file_type().is_symlink()
            || after.len() != self.expected.bytes
            || opened.len() != self.expected.bytes
            || self.bytes != self.expected.bytes
            || after.modified().ok() != self.modified
            || opened.modified().ok() != self.modified
            || same_file::Handle::from_path(&self.expected.path)
                .map_err(|e| fail(&e.to_string()))?
                != self.owner
            || hash != self.expected.sha256
        {
            return Err(fail("file_hash_or_change"));
        }
        Ok(())
    }
}
impl Read for Verified {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        let count = self.file.read(bytes)?;
        self.bytes = self
            .bytes
            .checked_add(count as u64)
            .filter(|n| *n <= self.expected.bytes)
            .ok_or_else(|| std::io::Error::other("observed_sam_intake_file_grew"))?;
        self.hash.update(&bytes[..count]);
        Ok(count)
    }
}

fn json_charge(bytes: u64) -> Result<usize, EngineError> {
    usize::try_from(bytes)
        .ok()
        .and_then(|n| n.checked_mul(64))
        .and_then(|n| n.checked_add(65536))
        .ok_or_else(|| fail("memory_budget"))
}
fn json<T: serde::de::DeserializeOwned>(
    file: &PrimerPayloadFile,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<T, EngineError> {
    let bytes = verified_bytes(file, file.bytes, check)?;
    serde_json::from_slice(&bytes).map_err(|e| fail(&format!("json: {e}")))
}

/// The same byte-bound reader for bounded locus FASTA and other explicit evidence inputs.
pub(crate) fn verified_bytes(
    file: &PrimerPayloadFile,
    cap: u64,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Vec<u8>, EngineError> {
    if file.bytes > cap {
        return Err(fail("input_byte_budget"));
    }
    let mut source = Verified::open(file)?;
    let mut bytes =
        Vec::with_capacity(usize::try_from(file.bytes).map_err(|_| fail("memory_budget"))?);
    let mut chunk = [0u8; 65536];
    loop {
        check()?;
        let n = source.read(&mut chunk).map_err(|e| fail(&e.to_string()))?;
        if n == 0 {
            break;
        }
        bytes.extend_from_slice(&chunk[..n]);
    }
    source.finish()?;
    Ok(bytes)
}

pub(crate) fn verified_json<T: serde::de::DeserializeOwned>(
    file: &PrimerPayloadFile,
    memory: usize,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<T, EngineError> {
    if json_charge(file.bytes)? > memory {
        return Err(fail("memory_budget"));
    }
    json(file, check)
}

pub(super) fn prepare(
    request: &Request,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Prepared, EngineError> {
    check()?;
    request
        .kernel_budgets
        .validate(request.required_flank_bases)?;
    let total = request
        .sam
        .bytes
        .checked_add(request.manifest.bytes)
        .and_then(|n| n.checked_add(request.dictionary.bytes))
        .ok_or_else(|| fail("input_budget"))?;
    if total > request.limits.max_input_bytes as u64
        || total > 128 * 1024 * 1024
        || !(1..=4096).contains(&request.max_units)
        || request.max_units > request.kernel_budgets.max_records
        || request.max_units > request.kernel_budgets.max_queries
        || request.limits.max_retained_bytes > 128 * 1024 * 1024
        || !text(&request.reference_id)
        || !text(&request.reference.assembly)
        || !sha(&request.reference.fasta_sha256)
    {
        return Err(fail("admission"));
    }
    for file in [&request.sam, &request.manifest, &request.dictionary] {
        path(&file.path)?;
        if !sha(&file.sha256) {
            return Err(fail("file_binding"));
        }
    }
    let valid_provenance = match &request.conversion {
        Conversion::DirectSam {
            producer,
            version,
            arguments,
            selection_scope,
        } => provenance(producer, version, arguments, selection_scope),
        Conversion::BamConversion {
            bam_sha256,
            producer,
            version,
            arguments,
            selection_scope,
        } => sha(bam_sha256) && provenance(producer, version, arguments, selection_scope),
        Conversion::ExecutedIndexedSelection { receipt } => receipt.validate(&request.sam).is_ok(),
    };
    if !valid_provenance {
        return Err(fail("conversion_provenance"));
    }
    let mut annotation_charge = 4096usize;
    let mut annotation_rows = 0usize;
    if request.annotation_sources.len() > 100_000 || request.comparison_genes.len() > 100_000 {
        return Err(fail("annotation_budget"));
    }
    for gene in std::iter::once(&request.reference_gene).chain(&request.comparison_genes) {
        check()?;
        if !identifier(&gene.source_id) || !identifier(&gene.gene_id) {
            return Err(fail("gene_selector"));
        }
        annotation_charge = annotation_charge
            .checked_add(16384)
            .ok_or_else(|| fail("memory_budget"))?;
    }
    for source in &request.annotation_sources {
        check()?;
        if !identifier(&source.source_id)
            || !sha(&source.annotation_sha256)
            || source.provider.as_ref().is_some_and(|s| !text(s))
            || [
                &source.annotation.schema,
                &source.annotation.reference_id,
                &source.annotation.assembly,
                &source.annotation.release,
                &source.annotation.reference_fasta_sha256,
                &source.annotation.coordinate_convention,
            ]
            .iter()
            .any(|s| !text(s))
        {
            return Err(fail("annotation_metadata"));
        }
        // Bound retained source metadata even for empty/unresolved sources, not only models.
        annotation_charge = annotation_charge
            .checked_add(131072)
            .ok_or_else(|| fail("memory_budget"))?;
        for transcript in &source.annotation.transcripts {
            check()?;
            annotation_rows = annotation_rows
                .checked_add(1)
                .filter(|n| *n <= request.kernel_budgets.max_annotation_rows && *n <= 100_000)
                .ok_or_else(|| fail("annotation_budget"))?;
            if !identifier(&transcript.transcript_id)
                || !identifier(&transcript.target_name)
                || transcript
                    .genomic_contig
                    .as_ref()
                    .is_some_and(|s| !identifier(s))
                || transcript.gene_ids.len() > 100_000
            {
                return Err(fail("annotation_metadata"));
            }
            let mut identifiers = 32768usize;
            for gene in &transcript.gene_ids {
                check()?;
                if !identifier(gene) {
                    return Err(fail("annotation_metadata"));
                }
                identifiers = identifiers
                    .checked_add(64 + 2 * gene.len())
                    .ok_or_else(|| fail("memory_budget"))?;
            }
            annotation_charge = transcript
                .genomic_exons
                .len()
                .checked_mul(128)
                .and_then(|n| n.checked_add(identifiers))
                .and_then(|n| n.checked_add(2048))
                .and_then(|n| n.checked_add(annotation_charge))
                .ok_or_else(|| fail("memory_budget"))?;
        }
    }
    let charged = json_charge(request.manifest.bytes)?
        .checked_add(json_charge(request.dictionary.bytes)?)
        .and_then(|n| n.checked_add(annotation_charge))
        .ok_or_else(|| fail("memory_budget"))?;
    let mut decoder_limits = request.limits.clone();
    decoder_limits.max_retained_bytes = request
        .limits
        .max_retained_bytes
        .checked_sub(charged)
        .filter(|n| *n > 0)
        .ok_or_else(|| fail("memory_budget"))?;
    let dictionary: ReferenceContigDictionary = json(&request.dictionary, check)?;
    if dictionary.schema != REFERENCE_CONTIG_DICTIONARY_SCHEMA
        || dictionary.reference_id != request.reference_id
        || dictionary.assembly != request.reference.assembly
        || dictionary.reference_fasta_sha256 != request.reference.fasta_sha256
        || dictionary.contigs.is_empty()
        || dictionary.contigs.len() > request.limits.max_dictionary_entries
    {
        return Err(fail("dictionary_binding"));
    }
    let mut contigs = BTreeMap::new();
    for entry in dictionary.contigs {
        check()?;
        if entry.length == 0
            || entry.length > i32::MAX as u64
            || contigs.insert(entry.name, entry.length as usize).is_some()
        {
            return Err(fail("dictionary_binding"));
        }
    }
    let manifest: Manifest = json(&request.manifest, check)?;
    if manifest.schema != MANIFEST_SCHEMA
        || manifest.sam_sha256 != request.sam.sha256
        || manifest.units.len() > request.max_units
    {
        return Err(fail("manifest_binding"));
    }
    let semantic = digest(&manifest)?;
    let retained_manifest = manifest.clone();
    let mut source = BufReader::new(Verified::open(&request.sam)?);
    let mut decoder = sam_text::Consumer::with_version_policy(
        &contigs,
        &decoder_limits,
        request.version_policy,
        check,
    )?;
    loop {
        check()?;
        let mut line = Vec::new();
        let n = source
            .by_ref()
            .take(request.limits.max_line_bytes as u64 + 1)
            .read_until(b'\n', &mut line)
            .map_err(|e| fail(&e.to_string()))?;
        if n == 0 {
            break;
        }
        decoder.push_line(&line, check)?;
    }
    source.into_inner().finish()?;
    let decoded = decoder.finish(check)?;
    if decoded.raw_sha256 != request.sam.sha256 || decoded.raw_bytes as u64 != request.sam.bytes {
        return Err(fail("sam_binding"));
    }
    let count = decoded.rows.len();
    let mut work = decoded.work;
    let mut rows = decoded.rows.into_iter().map(Some).collect::<Vec<_>>();
    let mut ids = BTreeSet::new();
    let mut units = Vec::new();
    let mut missing = Vec::new();
    for unit in manifest.units {
        check()?;
        if !identifier(&unit.unit_id)
            || !ids.insert(unit.unit_id.clone())
            || unit.rows.len() > request.kernel_budgets.max_segments_per_query
            || unit.rows.len() > 256
            || unit.ends.iter().any(|e| e.len() > 5)
            || unit.unclipped_axis_length == Some(0)
        {
            return Err(fail("unit_binding"));
        }
        let mut assigned = Vec::new();
        let mut axis_length = unit.unclipped_axis_length;
        for witness in unit.rows {
            check()?;
            let index = usize::try_from(witness.ordinal)
                .ok()
                .and_then(|n| n.checked_sub(1))
                .filter(|i| *i < count)
                .ok_or_else(|| fail("row_assignment"))?;
            let row = rows[index].take().ok_or_else(|| fail("row_assignment"))?;
            if row.row.ordinal != witness.ordinal
                || row.row.qname != witness.qname
                || row.read_group != witness.read_group
            {
                return Err(fail("row_witness"));
            }
            if assigned
                .first()
                .is_some_and(|previous: &sam_text::DecodedRow| {
                    previous.read_group != row.read_group
                })
                || (row.row.qname != "*"
                    && assigned
                        .iter()
                        .find(|previous| previous.row.qname != "*")
                        .is_some_and(|previous| previous.row.qname != row.row.qname))
            {
                return Err(fail("unit_row_identity"));
            }
            if let Some(observed) = sam::normalized_axis_length(&row.row)? {
                if axis_length.is_some_and(|declared| declared != observed) {
                    return Err(fail("axis_binding"));
                }
                axis_length = Some(observed);
            }
            assigned.push(row);
        }
        let mut complete_views = true;
        for row in &assigned {
            for companion in &row.companions {
                check()?;
                let mut found = false;
                for candidate in &assigned {
                    check()?;
                    work = work
                        .checked_add(companion.cigar.len())
                        .and_then(|n| n.checked_add(1))
                        .filter(|n| *n <= request.limits.max_work)
                        .ok_or_else(|| fail("work_budget"))?;
                    if candidate.row.ordinal != row.row.ordinal
                        && candidate.row.target.as_ref() == Some(&companion.target)
                        && candidate.row.position_1based == companion.position_1based
                        && candidate.row.mapq == companion.mapq
                        && (if candidate.row.flag & 16 == 0 {
                            '+'
                        } else {
                            '-'
                        }) == companion.strand
                        && candidate.row.cigar.as_ref() == Some(&companion.cigar)
                    {
                        found = true;
                        break;
                    }
                }
                if !found {
                    complete_views = false;
                    missing.push(MissingCompanion {
                        unit_id: unit.unit_id.clone(),
                        row_ordinal: row.row.ordinal,
                        companion: companion.clone(),
                    });
                }
            }
        }
        let masks = if let Some(witness) = unit.masks {
            if witness.coordinate_convention != "zero_based_half_open_unclipped_query"
                || Some(witness.unclipped_axis_length) != unit.unclipped_axis_length
                || !sha(&witness.source_sha256)
                || !text(&witness.derivation)
                || witness.intervals.len() > request.limits.max_operations
                || witness
                    .intervals
                    .iter()
                    .any(|&(a, z)| a >= z || z > witness.unclipped_axis_length)
                || !witness.intervals.windows(2).all(|w| w[0].1 < w[1].0)
            {
                return Err(fail("mask_witness"));
            }
            Some(witness.intervals)
        } else {
            None
        };
        let (search_complete, declared_views_complete) = if let Some(search) = unit.search {
            if !provenance(
                &search.producer,
                &search.version,
                &search.arguments,
                &search.selection_scope,
            ) || search.limitations.len() > 1024
                || search.limitations.iter().any(|s| !text(s))
            {
                return Err(fail("search_declaration"));
            }
            (
                search.search_complete && search.limitations.is_empty(),
                search.views_complete,
            )
        } else {
            (false, false)
        };
        units.push(sam::Unit {
            unit_id: unit.unit_id,
            declared_axis_length: axis_length,
            masks,
            ends: unit.ends,
            search_declared_complete: search_complete,
            views_declared_complete: declared_views_complete && complete_views,
            rows: assigned.into_iter().map(|r| r.row).collect(),
        });
    }
    if rows.iter().any(Option::is_some) {
        return Err(fail("row_unassigned"));
    }
    check()?;
    let receipt = Receipt {
        sam_sha256: request.sam.sha256.clone(),
        manifest_sha256: request.manifest.sha256.clone(),
        manifest_semantic_sha256: semantic,
        dictionary_sha256: request.dictionary.sha256.clone(),
        source_eof_verified: true,
        coordinate_authority: "producer_declared".into(),
        identity_scope: "sam_declared_query_unit".into(),
        version_policy: decoded.version_policy,
        unpadded_reference: request.unpadded_reference,
        declared_version: decoded.declared_version,
        counts: decoded.counts,
        query_units: units.len(),
        query_unit_manifest: retained_manifest,
        missing_sa_companions: missing,
        conversion: request.conversion.clone(),
        retained_peak_bytes: charged + decoded.retained_peak_bytes,
        work,
    };
    Ok(Prepared {
        input: sam::SamInput {
            source_sha256: request.sam.sha256.clone(),
            reference: request.reference.clone(),
            reference_gene: request.reference_gene.clone(),
            comparison_genes: request.comparison_genes.clone(),
            annotation_sources: request.annotation_sources.clone(),
            dictionary: contigs,
            unpadded_reference: request.unpadded_reference,
            units,
            required_flank_bases: request.required_flank_bases,
            budgets: request.kernel_budgets.clone(),
            max_rows: request.limits.max_rows,
        },
        receipt,
    })
}
