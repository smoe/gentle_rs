//! C1: authenticated coding-model derivation, without allele calls or HGVS interpretation.
//! Inputs are local immutable byte bindings; completeness is evidence, not provider inference.

use super::CdnaLocusReference;
use crate::engine::{EngineError, OperationProgress, PrimerPayloadAnnotation, PrimerPayloadFile};
pub use crate::observed_transcript_structures::annotation::VersionState as CodingVersionState;
use crate::observed_transcript_structures::{annotation::gtf, geometry, read_ends};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

mod model;
#[cfg(test)]
mod tests;

pub const CODING_MODELS_REQUEST_SCHEMA: &str = "gentle.coding_transcript_models_request.v1";
pub const CODING_MODELS_SCHEMA: &str = "gentle.coding_transcript_models.v1";

fn fail(message: impl std::fmt::Display) -> EngineError {
    EngineError::invalid_input(format!("coding_models.{message}"))
}
fn digest(bytes: &[u8]) -> String {
    crate::digest_utils::sha256_prefixed_bytes(bytes)
}
fn hash<T: Serialize>(value: &T) -> Result<String, EngineError> {
    Ok(digest(
        serde_json::to_value(value)
            .map_err(fail)?
            .to_string()
            .as_bytes(),
    ))
}
fn token(value: &str) -> bool {
    !value.is_empty() && value.len() <= 512 && !value.chars().any(char::is_control)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CdsRowConvention {
    StopCodonExcludedFromCdsRows,
    StopCodonIncludedInCdsRows,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CodingModelProvider {
    pub name: String,
    pub release: String,
    pub cds_row_convention: CdsRowConvention,
}

/// Explicit declaration; table 1 alone is not permission to interpret organellar sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CodingReferenceKind {
    Nuclear,
    Mitochondrial,
    Other,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CodingModelSelector {
    pub source_id: String,
    pub gene_id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CodingModelBudgets {
    pub max_input_bytes: u64,
    pub max_line_bytes: usize,
    pub max_rows: usize,
    pub max_transcripts: usize,
    pub max_segments: usize,
    pub max_transcript_bases: usize,
    pub max_work: usize,
    pub max_output_bytes: usize,
}
impl Default for CodingModelBudgets {
    fn default() -> Self {
        Self {
            max_input_bytes: 64_000_000,
            max_line_bytes: 65_536,
            max_rows: 500_000,
            max_transcripts: 256,
            max_segments: 100_000,
            max_transcript_bases: 2_000_000,
            max_work: 100_000_000,
            max_output_bytes: 8_000_000,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CodingTranscriptModelsRequest {
    pub schema: String,
    pub gtf: PrimerPayloadFile,
    pub normalized_annotation: PrimerPayloadFile,
    pub reference: CdnaLocusReference,
    pub selector: CodingModelSelector,
    pub provider: CodingModelProvider,
    pub reference_kind: CodingReferenceKind,
    pub translation_table: usize,
    #[serde(default)]
    pub protein_fasta: Option<PrimerPayloadFile>,
    #[serde(default)]
    pub budgets: CodingModelBudgets,
}

/// Accessions are withheld when versions are absent, conflicting or only partly declared.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CodingAccession {
    pub id: String,
    pub version_state: CodingVersionState,
    pub versioned_accession: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CodingSegment {
    pub genomic_0based_half_open: (usize, usize),
    pub gtf_frame: Option<u8>,
}

/// Geometry is available even when translation is withheld; normalized CDS includes a known stop.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CodingGeometry {
    pub segments: Vec<CodingSegment>,
    pub start_codon_segments: Vec<(usize, usize)>,
    pub stop_codon_segments: Vec<(usize, usize)>,
    pub utr_segments: BTreeMap<String, Vec<(usize, usize)>>,
    pub normalized_coding_cdna_0based_half_open: (usize, usize),
    pub normalized_genomic_coding_segments: Vec<(usize, usize)>,
    pub cdna_cds_start_1based: usize,
    pub cdna_cds_end_1based: usize,
    pub start_codon: Option<String>,
    pub stop_codon: Option<String>,
    pub cds_completeness_reasons: Vec<String>,
    pub protein: Option<CodingAccession>,
    pub translated_protein_sha256: Option<String>,
    pub supplied_protein_comparison: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CodingTranscriptModel {
    pub transcript: CodingAccession,
    pub gene_id: String,
    pub biotypes: Vec<String>,
    pub strand: char,
    pub genomic_exons_0based_half_open: Vec<(usize, usize)>,
    pub genomic_utr_segments: BTreeMap<String, Vec<(usize, usize)>>,
    pub transcript_completeness_tags: Vec<String>,
    pub tags: Vec<String>,
    pub coding: Option<CodingGeometry>,
    pub translation_status: String,
    pub translation_reasons: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CodingTranscriptModelsReport {
    pub schema: String,
    /// SHA-256 of canonical serde_json::Value text with this field omitted.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub content_sha256: String,
    pub request_sha256: String,
    pub request: CodingTranscriptModelsRequest,
    pub sequence_sha256: String,
    pub genetic_code_sha256: String,
    pub models: Vec<CodingTranscriptModel>,
    pub work: usize,
    pub limitations: Vec<String>,
}

struct Meter<'a> {
    remaining: usize,
    initial: usize,
    check: &'a mut dyn FnMut() -> Result<(), EngineError>,
}
impl Meter<'_> {
    fn charge(&mut self, amount: usize) -> Result<(), EngineError> {
        (self.check)()?;
        self.remaining = self
            .remaining
            .checked_sub(amount)
            .ok_or_else(|| fail("work_budget"))?;
        Ok(())
    }
}

#[derive(Default)]
struct Rows {
    strand: Option<char>,
    transcript_spans: BTreeSet<(usize, usize)>,
    exons: BTreeSet<(usize, usize)>,
    transcript_versions: BTreeSet<Option<String>>,
    protein_ids: BTreeSet<String>,
    protein_versions: BTreeSet<Option<String>>,
    biotypes: BTreeSet<String>,
    tags: BTreeSet<String>,
    cds: BTreeMap<(usize, usize), Option<u8>>,
    start: BTreeSet<(usize, usize)>,
    stop: BTreeSet<(usize, usize)>,
    utr: BTreeMap<String, BTreeSet<(usize, usize)>>,
}

fn bytes(
    file: &PrimerPayloadFile,
    limit: u64,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Vec<u8>, EngineError> {
    crate::observed_transcript_structures::indexed_locus::local_file(&file.path)?;
    crate::observed_transcript_structures::read_sam_bound_bytes(file, limit, check)
}

/// Return one authenticated, deterministic report. No project access, writes or provider calls.
pub fn derive_coding_transcript_models(
    request: CodingTranscriptModelsRequest,
    progress: &mut dyn FnMut(OperationProgress) -> bool,
) -> Result<CodingTranscriptModelsReport, EngineError> {
    let b = &request.budgets;
    let reference_hash = request
        .reference
        .genomic_reference_sha256
        .strip_prefix("sha256:");
    if request.schema != CODING_MODELS_REQUEST_SCHEMA
        || request.translation_table != 1
        || request.reference_kind != CodingReferenceKind::Nuclear
        || reference_hash.is_none_or(|h| h.len() != 64 || !h.bytes().all(|c| c.is_ascii_hexdigit()))
        || !token(&request.selector.source_id)
        || !token(&request.selector.gene_id)
        || !token(&request.reference.assembly)
        || !token(&request.reference.contig)
        || !token(&request.provider.name)
        || !token(&request.provider.release)
        || !(1..=128_000_000).contains(&b.max_input_bytes)
        || !(1..=1_048_576).contains(&b.max_line_bytes)
        || !(1..=1_000_000).contains(&b.max_rows)
        || !(1..=4096).contains(&b.max_transcripts)
        || !(1..=1_000_000).contains(&b.max_segments)
        || !(1..=16_000_000).contains(&b.max_transcript_bases)
        || !(1..=100_000_000).contains(&b.max_work)
        || !(1..=64_000_000).contains(&b.max_output_bytes)
    {
        return Err(fail("schema_context_or_budget"));
    }
    let total = [
        &request.gtf,
        &request.normalized_annotation,
        &request.reference.locus_fasta,
    ]
    .into_iter()
    .chain(request.protein_fasta.iter())
    .try_fold(0u64, |n, f| {
        n.checked_add(f.bytes).ok_or_else(|| fail("input_budget"))
    })?;
    if total > b.max_input_bytes {
        return Err(fail("input_budget"));
    }
    let mut check = || {
        if progress(OperationProgress::Workflow {
            completed: 0,
            total: 1,
        }) {
            Ok(())
        } else {
            Err(fail("cancelled"))
        }
    };
    let sequence = super::comparison::locus_sequence(
        &request.reference.locus_fasta,
        &request.reference.contig,
        &mut check,
    )?;
    let annotation: PrimerPayloadAnnotation = serde_json::from_slice(&bytes(
        &request.normalized_annotation,
        b.max_input_bytes,
        &mut check,
    )?)
    .map_err(fail)?;
    if annotation.schema != crate::engine::PRIMER_PAYLOAD_ANNOTATION_SCHEMA
        || annotation.reference_id != request.selector.source_id
        || annotation.assembly != request.reference.assembly
        || annotation.release != request.provider.release
        || annotation.reference_fasta_sha256 != request.reference.genomic_reference_sha256
        || annotation.coordinate_convention != "zero_based_half_open"
    {
        return Err(fail("annotation_binding"));
    }
    let raw = bytes(&request.gtf, b.max_input_bytes, &mut check)?;
    let protein_raw = request
        .protein_fasta
        .as_ref()
        .map(|f| bytes(f, b.max_input_bytes, &mut check))
        .transpose()?;
    let mut meter = Meter {
        remaining: b.max_work,
        initial: b.max_work,
        check: &mut check,
    };
    let mut selected = BTreeMap::new();
    for t in &annotation.transcripts {
        meter.charge(1)?;
        if t.gene_ids.contains(&request.selector.gene_id) {
            if t.genomic_exons.len() > b.max_segments
                || selected.insert(t.transcript_id.clone(), t).is_some()
                || selected.len() > b.max_transcripts
            {
                return Err(fail("transcript_universe"));
            }
        }
    }
    if selected.is_empty() {
        return Err(fail("transcript_universe"));
    }
    let mut rows = BTreeMap::<String, Rows>::new();
    let mut segments = 0usize;
    for (index, line) in std::str::from_utf8(&raw).map_err(fail)?.lines().enumerate() {
        meter.charge(line.len().saturating_add(1))?;
        if index >= b.max_rows || line.len() > b.max_line_bytes {
            return Err(fail("row_budget"));
        }
        if line.starts_with('#') {
            continue;
        }
        if line.chars().any(|c| c.is_control() && c != '\t') {
            return Err(fail("gtf_control"));
        }
        let row = gtf::parse_coding(line)?;
        if row.gene != request.selector.gene_id {
            if row
                .transcript
                .as_ref()
                .is_some_and(|id| selected.contains_key(id))
            {
                return Err(fail("transcript_gene_conflict"));
            }
            continue;
        }
        let Some(id) = row.transcript.as_ref() else {
            continue;
        };
        if !selected.contains_key(id)
            || row.contig != request.reference.contig
            || !matches!(row.strand, '+' | '-')
        {
            return Err(fail("gtf_transcript_universe"));
        }
        segments += 1;
        if segments > b.max_segments {
            return Err(fail("segment_budget"));
        }
        let value = rows.entry(id.clone()).or_default();
        if value
            .strand
            .replace(row.strand)
            .is_some_and(|s| s != row.strand)
        {
            return Err(fail("strand_conflict"));
        }
        value.transcript_versions.insert(row.transcript_version);
        for (key, values) in &row.attributes {
            match key.as_str() {
                "tag" => value.tags.extend(values.iter().cloned()),
                "transcript_biotype" | "transcript_type" => {
                    value.biotypes.extend(values.iter().cloned())
                }
                "protein_id" => value.protein_ids.extend(values.iter().cloned()),
                _ => (),
            }
        }
        if row.feature == "CDS" {
            value.protein_versions.insert(
                row.attributes
                    .get("protein_version")
                    .and_then(|v| v.first())
                    .cloned(),
            );
            if value.cds.insert(row.span, row.frame).is_some() {
                return Err(fail("duplicate_cds"));
            }
        } else {
            match row.feature {
                "transcript" => {
                    if !value.transcript_spans.insert(row.span) {
                        return Err(fail("duplicate_transcript"));
                    }
                }
                "exon" => {
                    if !value.exons.insert(row.span) {
                        return Err(fail("duplicate_exon"));
                    }
                }
                "start_codon" => {
                    if !value.start.insert(row.span) {
                        return Err(fail("duplicate_start"));
                    }
                }
                "stop_codon" => {
                    if !value.stop.insert(row.span) {
                        return Err(fail("duplicate_stop"));
                    }
                }
                "UTR" | "five_prime_utr" | "three_prime_utr" => {
                    value
                        .utr
                        .entry(row.feature.into())
                        .or_default()
                        .insert(row.span);
                }
                _ => (),
            }
        }
    }
    let proteins = model::proteins(protein_raw.as_deref(), b, &mut meter)?;
    let mut models = Vec::new();
    for (id, annotation) in selected {
        let row = rows
            .remove(&id)
            .ok_or_else(|| fail("missing_gtf_transcript"))?;
        models.push(model::derive(
            &request, &sequence, annotation, row, &proteins, &mut meter,
        )?);
    }
    let work = meter.initial - meter.remaining;
    let mut report = CodingTranscriptModelsReport {
        schema: CODING_MODELS_SCHEMA.into(), content_sha256: String::new(), request_sha256: hash(&request)?,
        request, sequence_sha256: digest(sequence.as_bytes()),
        genetic_code_sha256: digest(include_bytes!("../../assets/codon_tables.json")), models, work,
        limitations: vec!["Bound producer declarations are not provider authenticity or biological validation.".into(),
            "C1 derives coding models only; no HGVS, allele consequence, expression or clinical inference.".into()],
    };
    report.content_sha256 = hash(&report)?;
    if serde_json::to_vec(&report).map_err(fail)?.len() > report.request.budgets.max_output_bytes {
        return Err(fail("output_budget"));
    }
    Ok(report)
}
