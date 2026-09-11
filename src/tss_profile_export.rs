//! Report-only, staged TSS profile exports and their integrity verification.
//!
//! No sequence is read or rescored here. JSON retains the computational report;
//! TSV projects its window-start geometry and keeps raw and display values apart.
//! A flat, fresh staging directory is verified before publication. Receipts bind
//! every published file except themselves, and are integrity records, not signatures.

use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File, OpenOptions},
    io::{self, BufReader, BufWriter, Write},
    path::{Component, Path, PathBuf},
};

use gentle_engine::tss_profiles::{validate_panel, validate_record, validate_reference};
use gentle_protocol::{EngineError, ErrorCode, tss_profiles::*};
use gentle_render::tss_profiles::{TssRenderedPage, render_tss_profile_pages};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::{
    digest_utils::{sha256_file_hex, sha256_hex_bytes},
    svg_png::{SvgPngRenderOptions, SvgUsedFontIdentity},
};

const LOCKFILE: &str = include_str!("../Cargo.lock");
const INDEX_SCHEMA: &str = "gentle.tss_tfbs_profile_index.v1";
const REPORT_FILE: &str = "report.json";
const INDEX_FILE: &str = "index.json";
const RECEIPT_FILE: &str = "receipt.json";
const REQUEST_FILE: &str = "export-request.json";
const README_FILE: &str = "README.md";
const COMPARISON_FILE: &str = "comparisons.tsv";
const MAX_WINDOWS: usize = 4_096;
const MAX_MATRICES: usize = 256;
const MAX_MOTIF_BP: usize = 64;
const MAX_WINDOW_BP: usize = 1_000_000;
const MAX_SCORE_CELLS: usize = 10_000_000;
const MAX_TEXT_BYTES: usize = 16_384;
const MAX_BINDINGS: usize = 32_768;
const MAX_REPORT_BYTES: u64 = 256 * 1024 * 1024;
const MAX_METADATA_BYTES: u64 = 32 * 1024 * 1024;
const MAX_EXPORT_BYTES: u64 = 4 * 1024 * 1024 * 1024;
const MAX_OUTPUT_FILES: usize = 32_768;
const MAX_SVG_BYTES: usize = 32 * 1024 * 1024;
const MAX_RASTER_PIXELS: u64 = 40_000_000;
const FONT_IDENTITY_STATUS: &str = "glyph_used_font_sources_recorded";
const FONT_DIGEST_CONVENTION: &str =
    "sha256(complete font source/container bytes); face_index recorded separately";
const TARGET_BUNDLE_SCHEMA: &str = "gentle.target_tss_fasta_export.v1";
const GENERIC_SELECTION_LABEL: &str = "Selected TSS";
const GENERIC_SELECTION_LEGEND: &str = "Selected according to the supplied report. No descriptive evidence criterion was recorded; this does not establish TSS usage, direct binding or promoter activity.";

fn invalid(message: impl Into<String>) -> EngineError {
    EngineError::invalid_input(format!("TSS profile export: {}", message.into()))
}

fn io_error(context: &str, error: impl std::fmt::Display) -> EngineError {
    EngineError::new(
        ErrorCode::Io,
        format!("TSS profile export: {context}: {error}"),
    )
}

fn text_field(value: &str, field: &str) -> Result<(), EngineError> {
    if value.trim().is_empty() || value.len() > MAX_TEXT_BYTES || value.contains('\0') {
        return Err(invalid(format!(
            "{field} must be nonempty, bounded text without NUL"
        )));
    }
    Ok(())
}

fn digest_field(value: &str, field: &str) -> Result<(), EngineError> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(invalid(format!(
            "{field} must be a lowercase, unprefixed SHA-256"
        )));
    }
    Ok(())
}

fn portable_name(name: &str) -> bool {
    !name.is_empty()
        && name.len() <= 240
        && !name.starts_with('.')
        && !name.ends_with('.')
        && name
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b))
        && !matches!(
            name.split('.')
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
        )
}

fn validate_bindings(bindings: &[TssInputBinding]) -> Result<(), EngineError> {
    if bindings.len() > MAX_BINDINGS {
        return Err(invalid("too many input bindings"));
    }
    let mut seen = BTreeMap::new();
    for binding in bindings {
        text_field(&binding.role, "input role")?;
        text_field(&binding.name, "input name")?;
        let path = Path::new(&binding.name);
        if path.is_absolute()
            || binding.name.contains(['\\', ':'])
            || binding.name.chars().any(char::is_control)
            || path
                .components()
                .any(|c| !matches!(c, Component::Normal(_)))
        {
            return Err(invalid("input names must be portable relative paths"));
        }
        digest_field(&binding.sha256, "input digest")?;
        let key = (&binding.role, &binding.name);
        if let Some(previous) = seen.insert(key, &binding.sha256) {
            if previous != &binding.sha256 {
                return Err(invalid("conflicting input digests for one role/name"));
            }
        }
    }
    Ok(())
}

fn manifest_digest(bindings: &[TssInputBinding]) -> Result<&str, EngineError> {
    // The reader uses bundle_manifest; early portable reports used manifest.
    // Never infer this role from a filename or a substring in another role.
    let mut manifests = bindings
        .iter()
        .filter(|binding| matches!(binding.role.as_str(), "bundle_manifest" | "manifest"));
    let first = manifests
        .next()
        .ok_or_else(|| invalid("report/receipt requires an explicit input-manifest binding"))?;
    digest_field(&first.sha256, "input manifest digest")?;
    if manifests.any(|binding| binding.sha256 != first.sha256) {
        return Err(invalid("conflicting input-manifest digests"));
    }
    Ok(&first.sha256)
}

fn source_revision(report: &TssProfileReport) -> Option<&str> {
    report
        .source
        .as_ref()
        .and_then(|source| source.source_revision.as_deref())
}

fn provenance_text(value: &str, field: &str) -> Result<(), EngineError> {
    text_field(value, field)?;
    if value.trim() != value || value.chars().any(char::is_control) {
        return Err(invalid(format!(
            "{field} must be trimmed text without control characters"
        )));
    }
    Ok(())
}

fn validate_source(report: &TssProfileReport) -> Result<(), EngineError> {
    let manifest = manifest_digest(&report.inputs)?;
    if let Some(hash) = &report.producer_executable_sha256 {
        digest_field(hash, "profile producer executable digest")?;
    }
    if let Some(source) = &report.source {
        if !matches!(source.schema.as_str(), BUNDLE_SCHEMA | TARGET_BUNDLE_SCHEMA) {
            return Err(invalid("unsupported input source schema"));
        }
        digest_field(&source.manifest_sha256, "source manifest digest")?;
        if source.manifest_sha256 != manifest {
            return Err(invalid(
                "source manifest digest disagrees with the actual input-manifest binding",
            ));
        }
        for (name, value) in [
            ("source revision", &source.source_revision),
            ("source dataset ID", &source.dataset_id),
        ] {
            if let Some(value) = value {
                provenance_text(value, name)?;
            }
        }
        if let Some(hash) = &source.producer_sha256 {
            digest_field(hash, "bundle source producer digest")?;
        }
        if source.schema == TARGET_BUNDLE_SCHEMA
            && (source.source_revision.is_none() || source.dataset_id.is_none())
        {
            return Err(invalid(
                "target source requires its explicit source_revision and dataset_id",
            ));
        }
    }
    Ok(())
}

fn validate_selection_evidence(window: &TssProfileWindow) -> Result<(), EngineError> {
    if let Some(evidence) = &window.selection_evidence {
        if !window.selected {
            return Err(invalid(
                "selection evidence is only valid on a selected TSS",
            ));
        }
        for (name, value) in [
            ("selection label", &evidence.label),
            ("selected-panel legend", &evidence.legend),
            ("selection criterion", &evidence.criterion),
        ] {
            text_field(value, name)?;
        }
        if let Some(factor) = &evidence.factor {
            text_field(factor, "selection factor")?;
        }
    }
    Ok(())
}

fn selection_fields(window: &TssProfileWindow) -> [Option<&str>; 4] {
    match (&window.selection_evidence, window.selected) {
        (Some(evidence), true) => [
            Some(&evidence.label),
            Some(&evidence.legend),
            Some(&evidence.criterion),
            evidence.factor.as_deref(),
        ],
        (None, true) => [
            Some(GENERIC_SELECTION_LABEL),
            Some(GENERIC_SELECTION_LEGEND),
            None,
            None,
        ],
        _ => [None, None, None, None],
    }
}

struct SizeCounter {
    bytes: u64,
    limit: u64,
}

impl Write for SizeCounter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len() as u64)
            .filter(|n| *n <= self.limit)
            .ok_or_else(|| io::Error::other("serialized size limit exceeded"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

fn bounded_json(value: &impl Serialize, limit: u64) -> Result<(), EngineError> {
    serde_json::to_writer(SizeCounter { bytes: 0, limit }, value)
        .map_err(|e| invalid(format!("invalid or oversized JSON: {e}")))
}

fn validate_metadata(value: &Value) -> Result<(), EngineError> {
    let mut pending = vec![(value, 0)];
    let mut count = 0;
    while let Some((value, depth)) = pending.pop() {
        count += 1;
        if depth > 32 || count > 8_192 {
            return Err(invalid("normalization metadata exceeds depth/node limits"));
        }
        match value {
            Value::Array(values) => {
                if values.len() > 8_192 {
                    return Err(invalid("metadata array is too large"));
                }
                pending.extend(values.iter().map(|v| (v, depth + 1)));
            }
            Value::Object(values) => {
                if values.len() > 8_192 {
                    return Err(invalid("metadata object is too large"));
                }
                for (key, value) in values {
                    text_field(key, "metadata key")?;
                    pending.push((value, depth + 1));
                }
            }
            Value::String(text) if text.len() > MAX_TEXT_BYTES => {
                return Err(invalid("metadata string is too large"));
            }
            _ => {}
        }
    }
    bounded_json(value, 256 * 1024)
}

fn validate_geometry(geometry: &TssGeometry) -> Result<usize, EngineError> {
    gentle_engine::tss_profiles::validate_geometry(geometry)?;
    let length = geometry
        .length()
        .filter(|n| *n <= MAX_WINDOW_BP)
        .ok_or_else(|| invalid("window length exceeds its bound or overflows"))?;
    Ok(length)
}

fn validate_peak(peak: &TssPeak, scores: &[Option<f64>]) -> Result<(), EngineError> {
    if !peak.score.is_finite()
        || scores.get(peak.local_start_0based).copied().flatten() != Some(peak.score)
    {
        return Err(invalid(
            "peak must point to the same finite raw score in its track",
        ));
    }
    Ok(())
}

/// Recheck portable report invariants without scoring, rendering, or filesystem writes.
///
/// Score vectors contain exactly `max(0, window_length - motif_length + 1)`
/// entries per strand. A null entry is unavailable; trailing non-fitting starts
/// are not entries. Bounds apply to total *exported* cells, including those tails.
/// Every same-factor matrix pair requires both strand comparisons with counts
/// matching the supplied validity masks. Correlation values are not recomputed.
pub fn validate_tss_profile_report(report: &TssProfileReport) -> Result<(), EngineError> {
    if report.schema != REPORT_SCHEMA || report.panel_resolution.panel.schema != PANEL_SCHEMA {
        return Err(invalid("unsupported report or panel schema"));
    }
    if report.non_claims != NON_CLAIMS {
        return Err(invalid(
            "report non_claims must equal the canonical protocol non-claims",
        ));
    }
    if report.windows.is_empty() || report.windows.len() > MAX_WINDOWS {
        return Err(invalid("report must contain 1..=4096 TSS windows"));
    }
    validate_reference(&report.reference)?;
    for (name, value) in [
        ("producer_revision", &report.producer_revision),
        ("verification", &report.verification),
        ("non_claims", &report.non_claims),
    ] {
        text_field(value, name)?;
    }
    digest_field(&report.lockfile_sha256, "producer lockfile digest")?;
    if report.inputs.is_empty() || report.panel_resolution.registry_sources.is_empty() {
        return Err(invalid(
            "input and actual registry-source bindings are required",
        ));
    }
    validate_bindings(&report.inputs)?;
    validate_bindings(&report.panel_resolution.registry_sources)?;
    validate_source(report)?;
    if report.score_policy.is_empty()
        || report.score_policy.len() > 256
        || report.warnings.len() > MAX_WINDOWS
    {
        return Err(invalid("missing or oversized score policy/warnings"));
    }
    for (key, value) in &report.score_policy {
        text_field(key, "score policy key")?;
        text_field(value, "score policy value")?;
    }
    for warning in &report.warnings {
        text_field(warning, "warning")?;
    }
    let resolution = &report.panel_resolution;
    let panel = &resolution.panel;
    validate_panel(panel)?;
    digest_field(&resolution.panel_sha256, "panel digest")?;
    if let Some(url) = &resolution.registry_source_url {
        text_field(url, "registry source URL")?;
    }
    if panel.factors.len() > MAX_MATRICES || resolution.matrices.len() != panel.factors.len() {
        return Err(invalid("invalid panel, matrix count or top-hit limit"));
    }
    for (specification, matrix) in panel.factors.iter().zip(&resolution.matrices) {
        if specification.source_id.rsplit('.').next() != Some(matrix.version.as_str())
            || serde_json::to_value(specification).map_err(|e| invalid(e.to_string()))?
                != serde_json::to_value(&matrix.specification)
                    .map_err(|e| invalid(e.to_string()))?
        {
            return Err(invalid(
                "matrix identities, versions, ordering or score policies disagree",
            ));
        }
        let length = matrix.matrix_counts.len();
        if length == 0
            || length > MAX_MOTIF_BP
            || matrix.consensus.len() != length
            || !matrix.consensus.is_ascii()
        {
            return Err(invalid("invalid matrix/consensus shape"));
        }
        for column in &matrix.matrix_counts {
            let total = column.iter().sum::<f64>();
            if column.iter().any(|x| !x.is_finite() || *x < 0.0)
                || !total.is_finite()
                || total <= 0.0
            {
                return Err(invalid(
                    "matrix columns require finite nonnegative counts and a positive finite total",
                ));
            }
        }
        digest_field(&matrix.matrix_sha256, "matrix digest")?;
        let bytes = serde_json::to_vec(&(
            &specification.source_id,
            &Some(&specification.factor_id),
            &matrix.matrix_counts,
        ))
        .map_err(|e| invalid(e.to_string()))?;
        if sha256_hex_bytes(&bytes) != matrix.matrix_sha256 {
            return Err(invalid(
                "matrix digest does not match the actual accession, exact factor name and counts",
            ));
        }
    }
    let mut expected_pairs = BTreeSet::new();
    for (left, specification) in panel.factors.iter().enumerate() {
        for (right, other) in panel.factors.iter().enumerate().skip(left + 1) {
            if specification.factor_id == other.factor_id {
                for strand in [TssStrand::Plus, TssStrand::Minus] {
                    expected_pairs.insert((left, right, strand.as_str()));
                }
            }
        }
    }
    let mut promoters = BTreeSet::new();
    let mut genes = BTreeMap::new();
    let mut cells = 0usize;
    for window in &report.windows {
        let record = &window.record;
        validate_record(record)?;
        validate_selection_evidence(window)?;
        if !promoters.insert(&record.promoter_id) {
            return Err(invalid("duplicate promoter ID"));
        }
        if genes
            .insert(&record.gene_id, &record.gene_symbol)
            .is_some_and(|symbol| symbol != &record.gene_symbol)
        {
            return Err(invalid("one gene ID has conflicting symbols"));
        }
        let length = validate_geometry(&record.geometry)?;
        cells = length
            .checked_mul(panel.factors.len())
            .and_then(|n| n.checked_mul(2))
            .and_then(|n| cells.checked_add(n))
            .filter(|n| *n <= MAX_SCORE_CELLS)
            .ok_or_else(|| invalid("report exceeds the 10000000 exported score-cell limit"))?;
        if window.tracks.len() != resolution.matrices.len() {
            return Err(invalid("track count does not match the resolved panel"));
        }
        for (track, matrix) in window.tracks.iter().zip(&resolution.matrices) {
            let expected = length.saturating_sub(matrix.matrix_counts.len() - 1);
            if track.accession != matrix.specification.source_id
                || track.motif_length_bp != matrix.matrix_counts.len()
                || track.forward_scores.len() != expected
                || track.reverse_scores.len() != expected
            {
                return Err(invalid(
                    "track identity, motif length or score-vector shape disagrees with the report",
                ));
            }
            for (scores, maximum, peaks) in [
                (
                    &track.forward_scores,
                    &track.forward_maximum,
                    &track.forward_peaks,
                ),
                (
                    &track.reverse_scores,
                    &track.reverse_maximum,
                    &track.reverse_peaks,
                ),
            ] {
                if scores.iter().flatten().any(|x| !x.is_finite()) {
                    return Err(invalid("raw scores must be finite or null"));
                }
                if peaks.len() > panel.top_hit_count {
                    return Err(invalid("too many peaks"));
                }
                if let Some(peak) = maximum {
                    validate_peak(peak, scores)?;
                }
                let mut starts = BTreeSet::new();
                for peak in peaks {
                    validate_peak(peak, scores)?;
                    if !starts.insert(peak.local_start_0based) {
                        return Err(invalid("duplicate peak start"));
                    }
                }
            }
            validate_metadata(&track.normalization_reference)?;
        }
        if window.comparisons.len() != expected_pairs.len() {
            return Err(invalid(
                "comparisons must include every same-factor matrix pair on both local strands",
            ));
        }
        let mut pairs = BTreeSet::new();
        for comparison in &window.comparisons {
            text_field(&comparison.method, "comparison method")?;
            if let Some(reason) = &comparison.undefined_reason {
                text_field(reason, "undefined reason")?;
            }
            let left = panel
                .factors
                .iter()
                .position(|m| m.source_id == comparison.left_accession);
            let right = panel
                .factors
                .iter()
                .position(|m| m.source_id == comparison.right_accession);
            let (Some(left), Some(right)) = (left, right) else {
                return Err(invalid("comparison refers to an unknown matrix"));
            };
            let (a, b) = match comparison.local_strand {
                TssStrand::Plus => (
                    &window.tracks[left].forward_scores,
                    &window.tracks[right].forward_scores,
                ),
                TssStrand::Minus => (
                    &window.tracks[left].reverse_scores,
                    &window.tracks[right].reverse_scores,
                ),
            };
            // Count the supplied validity masks only. Neither scores nor
            // correlation coefficients are recomputed by an export.
            let paired_count = a
                .iter()
                .zip(b)
                .filter(|(a, b)| a.is_some() && b.is_some())
                .count();
            let excluded_count = a.len().max(b.len()) - paired_count;
            if left == right
                || panel.factors[left].factor_id != comparison.factor_id
                || panel.factors[right].factor_id != comparison.factor_id
                || !pairs.insert((
                    left.min(right),
                    left.max(right),
                    comparison.local_strand.as_str(),
                ))
                || comparison.paired_window_count != paired_count
                || comparison.excluded_window_count != excluded_count
                || [comparison.pearson, comparison.spearman]
                    .into_iter()
                    .flatten()
                    .any(|x| !x.is_finite() || !(-1.0..=1.0).contains(&x))
                || ((comparison.pearson.is_none() || comparison.spearman.is_none())
                    != comparison.undefined_reason.is_some())
                || ((comparison.pearson.is_some() || comparison.spearman.is_some())
                    && comparison.paired_window_count < 2)
            {
                return Err(invalid(
                    "invalid comparison identity, counts, outcome or finite correlation",
                ));
            }
        }
        if pairs != expected_pairs {
            return Err(invalid(
                "comparison pair/strand identities do not match the resolved panel",
            ));
        }
    }
    bounded_json(report, MAX_REPORT_BYTES)
}

fn validate_options(request: &ExportTssProfilesRequest) -> Result<(), EngineError> {
    if request.formats.is_empty()
        || request.formats.len() > 3
        || request
            .formats
            .iter()
            .enumerate()
            .any(|(i, f)| request.formats[..i].contains(f))
        || request.rendering.panels_per_page == 0
        || request.rendering.panels_per_page > 32
    {
        return Err(invalid(
            "choose unique SVG/PNG/PDF formats and 1..=32 panels per page",
        ));
    }
    Ok(())
}

fn validate_rendering(
    report: &TssProfileReport,
    request: &ExportTssProfilesRequest,
) -> Result<(), EngineError> {
    validate_options(request)?;
    let panel = &report.panel_resolution.panel;
    if request.rendering.scale_mode.unwrap_or(panel.scale_mode) == TssScaleMode::Shared
        && (panel.calibration_state != TssCalibrationState::CrossSourceCalibrated
            || panel.calibration_id.is_none()
            || panel.calibration_sha256.is_none())
    {
        return Err(invalid(
            "shared scales require a typed cross-source calibration ID and digest",
        ));
    }
    Ok(())
}

/// Check destination freshness and request options before expensive computation.
/// Does not create directories or files; the destination's parent must exist.
/// Calibration compatibility is checked later against the computed report.
pub fn preflight_tss_export(request: &ExportTssProfilesRequest) -> Result<(), EngineError> {
    validate_options(request)?;
    destination(&request.output_dir)?;
    Ok(())
}

fn checkpoint(should_continue: &mut dyn FnMut() -> bool) -> Result<(), EngineError> {
    if !should_continue() {
        return Err(invalid(
            "cancelled before publication; no success receipt was published",
        ));
    }
    Ok(())
}

struct Gene<'a> {
    id: &'a str,
    symbol: &'a str,
    stem: String,
    windows: Vec<&'a TssProfileWindow>,
}

fn gene_stem(symbol: &str, id: &str) -> String {
    fn slug(text: &str, limit: usize) -> String {
        let value: String = text
            .chars()
            .map(|c| {
                if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                    c
                } else {
                    '_'
                }
            })
            .take(limit)
            .collect();
        if value.is_empty() {
            "unnamed".into()
        } else {
            value
        }
    }
    let identity = format!("{}:{symbol}{}:{id}", symbol.len(), id.len());
    format!(
        "gene-{}--{}--{}",
        slug(symbol, 32),
        slug(id, 48),
        sha256_hex_bytes(identity.as_bytes())
    )
}

fn genes(report: &TssProfileReport) -> Result<Vec<Gene<'_>>, EngineError> {
    let mut genes: Vec<Gene<'_>> = Vec::new();
    let mut positions = BTreeMap::new();
    let mut names = BTreeSet::new();
    for window in &report.windows {
        let record = &window.record;
        let position = if let Some(position) = positions.get(record.gene_id.as_str()) {
            *position
        } else {
            let stem = gene_stem(&record.gene_symbol, &record.gene_id);
            if !portable_name(&stem) || !names.insert(stem.to_ascii_lowercase()) {
                return Err(invalid("gene output names collide or are unsafe"));
            }
            let position = genes.len();
            genes.push(Gene {
                id: &record.gene_id,
                symbol: &record.gene_symbol,
                stem,
                windows: Vec::new(),
            });
            positions.insert(record.gene_id.as_str(), position);
            position
        };
        genes[position].windows.push(window);
    }
    Ok(genes)
}

// Borrow the unchanged protocol fields; do not clone all score arrays to subset a gene.
#[derive(Serialize)]
struct GeneReport<'a> {
    schema: &'a str,
    reference: &'a TssReference,
    panel_resolution: &'a TssPanelResolution,
    inputs: &'a [TssInputBinding],
    windows: &'a [&'a TssProfileWindow],
    #[serde(skip_serializing_if = "Option::is_none")]
    source: Option<&'a TssBundleSource>,
    producer_revision: &'a str,
    #[serde(skip_serializing_if = "Option::is_none")]
    producer_executable_sha256: Option<&'a str>,
    lockfile_sha256: &'a str,
    score_policy: &'a BTreeMap<String, String>,
    verification: &'a str,
    warnings: &'a [String],
    non_claims: &'a str,
}

impl<'a> GeneReport<'a> {
    fn new(report: &'a TssProfileReport, gene: &'a Gene<'a>) -> Self {
        Self {
            schema: &report.schema,
            reference: &report.reference,
            panel_resolution: &report.panel_resolution,
            inputs: &report.inputs,
            windows: &gene.windows,
            source: report.source.as_ref(),
            producer_revision: &report.producer_revision,
            producer_executable_sha256: report.producer_executable_sha256.as_deref(),
            lockfile_sha256: &report.lockfile_sha256,
            score_policy: &report.score_policy,
            verification: &report.verification,
            warnings: &report.warnings,
            non_claims: &report.non_claims,
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PageIndex {
    page_number: usize,
    page_count: usize,
    promoter_ids: Vec<String>,
    files: Vec<String>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct GeneIndex {
    gene_id: String,
    gene_symbol: String,
    promoter_ids: Vec<String>,
    selections: Vec<SelectionIndex>,
    report: String,
    scores: String,
    comparisons: String,
    pages: Vec<PageIndex>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct SelectionIndex {
    promoter_id: String,
    label: String,
    legend: String,
    /// None remains not assessed, rather than inferring a CUT&RUN criterion.
    criterion: Option<String>,
    factor: Option<String>,
}

fn selection_index(gene: &Gene<'_>) -> Vec<SelectionIndex> {
    gene.windows
        .iter()
        .filter(|window| window.selected)
        .map(|window| {
            let [label, legend, criterion, factor] = selection_fields(window);
            SelectionIndex {
                promoter_id: window.record.promoter_id.clone(),
                label: label.unwrap_or(GENERIC_SELECTION_LABEL).into(),
                legend: legend.unwrap_or(GENERIC_SELECTION_LEGEND).into(),
                criterion: criterion.map(str::to_owned),
                factor: factor.map(str::to_owned),
            }
        })
        .collect()
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ExportIndex {
    schema: String,
    report_sha256: String,
    input_manifest_sha256: String,
    source_revision: Option<String>,
    source: Option<TssBundleSource>,
    reference: TssReference,
    producer_revision: String,
    producer_executable_sha256: Option<String>,
    producer_lockfile_sha256: String,
    score_policy: BTreeMap<String, String>,
    rendering: TssProfileRenderOptions,
    formats: Vec<TssExportFormat>,
    tss_count: usize,
    page_count: usize,
    genes: Vec<GeneIndex>,
    /// Data files only. The receipt binds this index; neither hashes the receipt.
    outputs: BTreeMap<String, String>,
}

fn tsv(text: &str) -> String {
    let mut escaped = String::with_capacity(text.len());
    if text.starts_with('#') {
        escaped.push('\\');
    }
    for c in text.chars() {
        match c {
            '\\' => escaped.push_str("\\\\"),
            '\t' => escaped.push_str("\\t"),
            '\r' => escaped.push_str("\\r"),
            '\n' => escaped.push_str("\\n"),
            _ => escaped.push(c),
        }
    }
    escaped
}

fn nullable(value: Option<impl std::fmt::Display>) -> String {
    value
        .map(|v| v.to_string())
        .unwrap_or_else(|| "null".into())
}

fn genomic_strand(transcript: TssStrand, local: TssStrand) -> TssStrand {
    if local == TssStrand::Plus {
        transcript
    } else {
        transcript.opposite()
    }
}

fn write_tsv_preamble(
    writer: &mut dyn Write,
    report: &TssProfileReport,
) -> Result<(), EngineError> {
    let policy = serde_json::to_string(&report.score_policy).map_err(|e| invalid(e.to_string()))?;
    writeln!(writer, "# non_claims: {NON_CLAIMS}\n\
# coordinates: transcript-oriented window starts, local 0-based; TSS-relative and genomic 1-based coordinates use TssGeometry.relative_at/genomic_at. Genomic motif intervals are inclusive and ascending. Local reverse motifs have a distinct 5-prime endpoint, not a different common-axis start.\n\
# unavailable: literal null is not zero; non-fitting motif starts have null interval/endpoint and raw/display values.\n\
# methods: export copies raw report values and comparisons without rescoring. Comparison method, strand, sample/exclusion counts and undefined reason are retained per row; display clipping never changes comparison input.\n\
# score_kind: {}\n\
# display_clip_negative: {} (display_score only; raw_score is unchanged)\n\
# score_policy_json: {policy}", report.panel_resolution.panel.score_kind, report.panel_resolution.panel.clip_negative)
        .map_err(|e| io_error("write TSV non-claims and methods preamble", e))?;
    writeln!(writer, "# input_manifest_sha256: {}\n# source_revision: {}\n# source_schema: {}\n# dataset_id: {}\n# bundle_source_producer_sha256: {}\n# annotation_release: {}\n# producer_revision: {}\n# producer_executable_sha256: {}\n# missing_provenance: null means not assessed; source, profile producer and exporter identities are distinct.\n# selection_non_claim: Selection and recorded descriptive criteria do not establish TSS usage, direct binding or promoter activity.",
        manifest_digest(&report.inputs)?, nullable(source_revision(report).map(tsv)),
        nullable(report.source.as_ref().map(|source| tsv(&source.schema))),
        nullable(report.source.as_ref().and_then(|source| source.dataset_id.as_deref()).map(tsv)),
        nullable(report.source.as_ref().and_then(|source| source.producer_sha256.as_deref())),
        nullable(report.reference.annotation_release.as_deref().map(tsv)), tsv(&report.producer_revision),
        nullable(report.producer_executable_sha256.as_deref()))
        .map_err(|e| io_error("write TSV source provenance", e))
}

fn selection_tsv(window: &TssProfileWindow) -> String {
    selection_fields(window)
        .map(|field| field.map(tsv).unwrap_or_else(|| "null".into()))
        .join("\t")
}

fn write_selection_preamble(
    writer: &mut dyn Write,
    windows: &[&TssProfileWindow],
) -> Result<(), EngineError> {
    for window in windows.iter().filter(|window| window.selected) {
        let [label, legend, criterion, factor] =
            selection_fields(window).map(|field| field.map(tsv).unwrap_or_else(|| "null".into()));
        writeln!(writer, "# selection_promoter_id: {}\n# selection_label: {label}\n# selection_legend: {legend}\n# selection_criterion: {criterion}\n# selection_factor: {factor}", tsv(&window.record.promoter_id))
            .map_err(|e| io_error("write selected-panel TSV legend", e))?;
    }
    Ok(())
}

fn write_scores(
    writer: &mut dyn Write,
    report: &TssProfileReport,
    windows: &[&TssProfileWindow],
    should_continue: &mut dyn FnMut() -> bool,
) -> Result<(), EngineError> {
    write_tsv_preamble(writer, report)?;
    write_selection_preamble(writer, windows)?;
    writeln!(writer, "gene_id\tgene_symbol\tpromoter_id\tchromosome\ttranscript_strand\tselected\taccession\tfactor_id\tscore_kind\tlocal_window_start_0based\ttss_relative_window_start_bp\tgenomic_window_start_1based\tmotif_length_bp\tmotif_genomic_start_1based\tmotif_genomic_end_1based\tmotif_5prime_genomic_1based\tlocal_motif_strand\tgenomic_motif_strand\tavailability\traw_score\tdisplay_score\tdisplay_clip_negative\tselection_label\tselection_legend\tselection_criterion\tselection_factor")
        .map_err(|e| io_error("write score TSV header", e))?;
    let panel = &report.panel_resolution.panel;
    for window in windows {
        checkpoint(should_continue)?;
        let record = &window.record;
        let geometry = &record.geometry;
        let selection = selection_tsv(window);
        let prefix = format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            tsv(&record.gene_id),
            tsv(&record.gene_symbol),
            tsv(&record.promoter_id),
            tsv(&geometry.chromosome),
            geometry.strand.as_str(),
            window.selected
        );
        let length = geometry
            .length()
            .ok_or_else(|| invalid("window length overflow"))?;
        for (track, specification) in window.tracks.iter().zip(&panel.factors) {
            let prefix = format!(
                "{prefix}\t{}\t{}\t{}",
                tsv(&track.accession),
                tsv(&specification.factor_id),
                tsv(&panel.score_kind)
            );
            for start in 0..length {
                if start % 4_096 == 0 {
                    checkpoint(should_continue)?;
                }
                let relative = geometry
                    .relative_at(start)
                    .ok_or_else(|| invalid("relative coordinate overflow"))?;
                let genomic = geometry
                    .genomic_at(start)
                    .ok_or_else(|| invalid("genomic coordinate overflow"))?;
                let last = start
                    .checked_add(track.motif_length_bp - 1)
                    .filter(|end| *end < length);
                let interval = last
                    .and_then(|end| geometry.genomic_at(end))
                    .map(|end| (genomic.min(end), genomic.max(end), end));
                for (local, scores) in [
                    (TssStrand::Plus, &track.forward_scores),
                    (TssStrand::Minus, &track.reverse_scores),
                ] {
                    let raw = scores.get(start).copied().flatten();
                    let display = raw.map(|value| {
                        if panel.clip_negative {
                            value.max(0.0)
                        } else {
                            value
                        }
                    });
                    let endpoint = interval.map(|(_, _, end)| {
                        if local == TssStrand::Plus {
                            genomic
                        } else {
                            end
                        }
                    });
                    let availability = if interval.is_none() {
                        "outside_window"
                    } else if raw.is_none() {
                        "unavailable"
                    } else {
                        "available"
                    };
                    writeln!(writer, "{prefix}\t{start}\t{relative}\t{genomic}\t{}\t{}\t{}\t{}\t{}\t{}\t{availability}\t{}\t{}\t{}\t{selection}",
                        track.motif_length_bp, nullable(interval.map(|v| v.0)), nullable(interval.map(|v| v.1)),
                        nullable(endpoint), local.as_str(), genomic_strand(geometry.strand, local).as_str(),
                        nullable(raw), nullable(display), panel.clip_negative)
                        .map_err(|e| io_error("stream score TSV", e))?;
                }
            }
        }
    }
    Ok(())
}

fn write_comparisons<'a>(
    writer: &mut dyn Write,
    report: &TssProfileReport,
    windows: impl IntoIterator<Item = &'a TssProfileWindow>,
) -> Result<(), EngineError> {
    let windows = windows.into_iter().collect::<Vec<_>>();
    write_tsv_preamble(writer, report)?;
    write_selection_preamble(writer, &windows)?;
    writeln!(writer, "gene_id\tgene_symbol\tpromoter_id\tfactor_id\tleft_accession\tright_accession\tlocal_strand\tgenomic_strand\tinput_score_kind\tpaired_window_count\texcluded_window_count\tpearson\tspearman\tundefined_reason\tmethod\tinput_values\tdisplay_clip_negative\tselected\tselection_label\tselection_legend\tselection_criterion\tselection_factor")
        .map_err(|e| io_error("write comparison TSV header", e))?;
    let panel = &report.panel_resolution.panel;
    for window in windows {
        let selection = selection_tsv(window);
        for comparison in &window.comparisons {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\traw_unclipped\t{}\t{}\t{selection}",
                tsv(&window.record.gene_id),
                tsv(&window.record.gene_symbol),
                tsv(&window.record.promoter_id),
                tsv(&comparison.factor_id),
                tsv(&comparison.left_accession),
                tsv(&comparison.right_accession),
                comparison.local_strand.as_str(),
                genomic_strand(window.record.geometry.strand, comparison.local_strand).as_str(),
                tsv(&panel.score_kind),
                comparison.paired_window_count,
                comparison.excluded_window_count,
                nullable(comparison.pearson),
                nullable(comparison.spearman),
                comparison
                    .undefined_reason
                    .as_deref()
                    .map(tsv)
                    .unwrap_or_else(|| "null".into()),
                tsv(&comparison.method),
                panel.clip_negative,
                window.selected
            )
            .map_err(|e| io_error("stream comparison TSV", e))?;
        }
    }
    Ok(())
}

/// Resolve only existing, ordinary directory components. In particular, do not
/// canonicalize away a caller's symlink and then mistakenly call it safe.
fn checked_directory(path: &Path) -> Result<PathBuf, EngineError> {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()
            .map_err(|e| io_error("read working directory", e))?
            .join(path)
    };
    let mut resolved = PathBuf::new();
    for component in absolute.components() {
        match component {
            Component::ParentDir => {
                return Err(invalid("parent traversal is not allowed in output paths"));
            }
            Component::CurDir => continue,
            _ => resolved.push(component.as_os_str()),
        }
        let metadata =
            fs::symlink_metadata(&resolved).map_err(|e| io_error("inspect output directory", e))?;
        if metadata.file_type().is_symlink() || !metadata.is_dir() {
            return Err(invalid(
                "output directory and its ancestors must be ordinary directories, not symlinks",
            ));
        }
    }
    Ok(resolved)
}

fn destination(path: &str) -> Result<PathBuf, EngineError> {
    if path.is_empty() || path.chars().any(char::is_control) {
        return Err(invalid(
            "output_dir must be a nonempty path without control characters",
        ));
    }
    let path = Path::new(path);
    if path.components().any(|c| matches!(c, Component::ParentDir)) {
        return Err(invalid("output_dir cannot contain parent traversal"));
    }
    let name = path
        .file_name()
        .and_then(|s| s.to_str())
        .filter(|s| portable_name(s))
        .ok_or_else(|| invalid("output directory name is not a safe portable filename"))?;
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let output = checked_directory(parent)?.join(name);
    check_fresh_destination(&output)?;
    Ok(output)
}

fn check_fresh_destination(path: &Path) -> Result<(), EngineError> {
    match fs::symlink_metadata(path) {
        Ok(metadata) => {
            if metadata.file_type().is_symlink() || !metadata.is_dir() {
                return Err(invalid("destination is a symlink or is not a directory"));
            }
            if fs::read_dir(path)
                .map_err(|e| io_error("inspect destination contents", e))?
                .next()
                .transpose()
                .map_err(|e| io_error("inspect destination entry", e))?
                .is_some()
            {
                return Err(invalid(
                    "destination is nonempty; export history is never overwritten",
                ));
            }
            Ok(())
        }
        Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(io_error("inspect destination", error)),
    }
}

struct BudgetWriter<W> {
    inner: W,
    counter: SizeCounter,
}

impl<W: Write> Write for BudgetWriter<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.counter.write(bytes)?;
        self.inner.write_all(bytes)?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

struct Inventory<'a> {
    directory: &'a Path,
    hashes: BTreeMap<String, String>,
    bytes: u64,
}

impl Inventory<'_> {
    fn reserve(&self, name: &str) -> Result<(PathBuf, File), EngineError> {
        if !portable_name(name)
            || self.hashes.contains_key(name)
            || self.hashes.len() >= MAX_OUTPUT_FILES
        {
            return Err(invalid("unsafe, duplicate, or excessive output filenames"));
        }
        let path = self.directory.join(name);
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .map_err(|e| io_error("create new staged artifact", e))?;
        Ok((path, file))
    }

    fn record(&mut self, name: &str, path: &Path) -> Result<(), EngineError> {
        let metadata =
            fs::symlink_metadata(path).map_err(|e| io_error("inspect staged artifact", e))?;
        if !metadata.is_file() || metadata.file_type().is_symlink() {
            return Err(invalid("staged artifact is not an ordinary file"));
        }
        self.bytes = self
            .bytes
            .checked_add(metadata.len())
            .filter(|n| *n <= MAX_EXPORT_BYTES)
            .ok_or_else(|| invalid("export exceeds its total byte limit"))?;
        self.hashes.insert(
            name.into(),
            sha256_file_hex(path).map_err(|e| io_error("hash staged artifact", e))?,
        );
        Ok(())
    }

    fn write(
        &mut self,
        name: &str,
        write: impl FnOnce(&mut dyn Write) -> Result<(), EngineError>,
    ) -> Result<(), EngineError> {
        let (path, file) = self.reserve(name)?;
        let mut writer = BudgetWriter {
            inner: BufWriter::new(file),
            counter: SizeCounter {
                bytes: 0,
                limit: MAX_EXPORT_BYTES - self.bytes,
            },
        };
        write(&mut writer)?;
        writer
            .flush()
            .map_err(|e| io_error("flush staged artifact", e))?;
        writer
            .inner
            .get_ref()
            .sync_all()
            .map_err(|e| io_error("sync staged artifact", e))?;
        drop(writer);
        self.record(name, &path)
    }

    fn json(&mut self, name: &str, value: &impl Serialize) -> Result<(), EngineError> {
        self.write(name, |writer| {
            serde_json::to_writer(writer, value).map_err(|e| io_error("serialize artifact", e))
        })
    }

    fn bytes(&mut self, name: &str, bytes: &[u8]) -> Result<(), EngineError> {
        self.write(name, |writer| {
            writer
                .write_all(bytes)
                .map_err(|e| io_error("write artifact", e))
        })
    }
}

fn input_bindings(report: &TssProfileReport) -> Result<Vec<TssInputBinding>, EngineError> {
    let mut bindings = BTreeMap::new();
    for binding in report
        .inputs
        .iter()
        .chain(&report.panel_resolution.registry_sources)
        .chain(
            report
                .windows
                .iter()
                .filter_map(|w| w.detail_context.as_ref())
                .flat_map(|c| &c.bindings),
        )
    {
        let key = (binding.role.clone(), binding.name.clone());
        if bindings
            .insert(key, binding.sha256.clone())
            .is_some_and(|hash| hash != binding.sha256)
        {
            return Err(invalid("conflicting report and registry input bindings"));
        }
    }
    // Matrices bind [accession, exact factor name, counts]; normalized sequences
    // have their own identity, distinct from the exact source-file bytes.
    let mut add = |role: &str, name: String, hash: &str| -> Result<(), EngineError> {
        if bindings
            .insert((role.into(), name), hash.into())
            .is_some_and(|previous| previous != hash)
        {
            return Err(invalid("conflicting derived input binding"));
        }
        Ok(())
    };
    add(
        "resolved_panel",
        "panel.json".into(),
        &report.panel_resolution.panel_sha256,
    )?;
    for matrix in &report.panel_resolution.matrices {
        add(
            "resolved_matrix",
            format!("{}.json", matrix.specification.source_id),
            &matrix.matrix_sha256,
        )?;
    }
    for window in &report.windows {
        add(
            "normalized_sequence",
            format!(
                "sequence-{}.txt",
                sha256_hex_bytes(window.record.promoter_id.as_bytes())
            ),
            &window.record.sequence_sha256,
        )?;
    }
    if let Some(hash) = &report.panel_resolution.panel.calibration_sha256 {
        add("cross_matrix_calibration", "calibration.json".into(), hash)?;
    }
    if let Some(hash) = &report.producer_executable_sha256 {
        add(
            "producer_executable",
            "profile-producer-executable".into(),
            hash,
        )?;
    }
    if let Some(hash) = report
        .source
        .as_ref()
        .and_then(|source| source.producer_sha256.as_deref())
    {
        add(
            "bundle_source_producer",
            "bundle-source-producer".into(),
            hash,
        )?;
    }
    let result = bindings
        .into_iter()
        .map(|((role, name), sha256)| TssInputBinding { role, name, sha256 })
        .collect::<Vec<_>>();
    validate_bindings(&result)?;
    manifest_digest(&result)?;
    Ok(result)
}

fn svg_dimensions(svg: &str) -> Result<(u32, u32), EngineError> {
    let mut reader = quick_xml::Reader::from_str(svg);
    loop {
        match reader
            .read_event()
            .map_err(|e| invalid(format!("invalid rendered SVG: {e}")))?
        {
            quick_xml::events::Event::Start(root) | quick_xml::events::Event::Empty(root) => {
                if root.name().as_ref() != "svg" {
                    return Err(invalid("renderer did not return an SVG root"));
                }
                let mut width = None;
                let mut height = None;
                for attribute in root.attributes() {
                    let attribute = attribute.map_err(|e| invalid(e.to_string()))?;
                    let target = match attribute.key.as_ref() {
                        "width" => &mut width,
                        "height" => &mut height,
                        _ => continue,
                    };
                    let text = attribute.value.as_ref();
                    let value: f64 = text
                        .trim()
                        .strip_suffix("px")
                        .unwrap_or(text.trim())
                        .parse()
                        .map_err(|_| invalid("SVG dimensions must be explicit CSS pixels"))?;
                    if !value.is_finite() || value <= 0.0 || value > 65_535.0 {
                        return Err(invalid("rendered SVG dimensions exceed the export budget"));
                    }
                    *target = Some(value.ceil() as u32);
                }
                return width
                    .zip(height)
                    .ok_or_else(|| invalid("SVG lacks explicit bounded width/height"));
            }
            quick_xml::events::Event::Decl(_) | quick_xml::events::Event::Comment(_) => {}
            quick_xml::events::Event::Text(text)
                if text.as_ref().bytes().all(|b| b.is_ascii_whitespace()) => {}
            _ => return Err(invalid("unexpected content before SVG root")),
        }
    }
}

fn validate_pages(
    pages: &[TssRenderedPage],
    genes: &[Gene<'_>],
    raster: bool,
) -> Result<(), EngineError> {
    if pages.is_empty() || pages.len() > MAX_OUTPUT_FILES / 3 {
        return Err(invalid("renderer returned an empty or excessive page set"));
    }
    let mut total = 0usize;
    let mut seen_pages = BTreeSet::new();
    let mut seen_promoters = BTreeSet::new();
    let mut counts = BTreeMap::new();
    for page in pages {
        let gene = genes
            .iter()
            .find(|g| g.id == page.gene_id)
            .ok_or_else(|| invalid("page belongs to an unknown gene"))?;
        if gene.symbol != page.gene_symbol
            || page.page_number == 0
            || page.page_number > page.page_count
            || !seen_pages.insert((gene.id, page.page_number))
            || page.promoter_ids.is_empty()
            || page.svg.len() > MAX_SVG_BYTES
        {
            return Err(invalid(
                "invalid rendered page identity, numbering, membership or size",
            ));
        }
        let mut members = BTreeSet::new();
        for id in &page.promoter_ids {
            if !members.insert(id) || !gene.windows.iter().any(|w| w.record.promoter_id == *id) {
                return Err(invalid("page promoter membership disagrees with its gene"));
            }
            seen_promoters.insert(id.as_str());
        }
        let (width, height) = svg_dimensions(&page.svg)?;
        if raster && u64::from(width) * u64::from(height) > MAX_RASTER_PIXELS {
            return Err(invalid("page exceeds the 40000000-pixel raster budget"));
        }
        total = total
            .checked_add(page.svg.len())
            .filter(|n| *n <= MAX_REPORT_BYTES as usize)
            .ok_or_else(|| invalid("total rendered SVG size exceeds the export budget"))?;
        *counts.entry(gene.id).or_insert(0usize) += 1;
    }
    for page in pages {
        if counts.get(page.gene_id.as_str()).copied() != Some(page.page_count) {
            return Err(invalid("renderer page counts are inconsistent"));
        }
    }
    if genes
        .iter()
        .flat_map(|g| &g.windows)
        .any(|w| !seen_promoters.contains(w.record.promoter_id.as_str()))
    {
        return Err(invalid("renderer omitted a TSS window"));
    }
    Ok(())
}

// Cargo.lock may contain several resvg versions. Read the root's dependency
// reference rather than reporting whichever package happens to come first.
fn root_dependency_version(dependency: &str) -> Option<String> {
    let blocks = LOCKFILE.split("[[package]]").collect::<Vec<_>>();
    let root = blocks
        .iter()
        .find(|block| block.lines().any(|line| line.trim() == "name = \"GENtle\""))?;
    let reference = root
        .lines()
        .filter_map(|line| line.trim().strip_prefix('"')?.strip_suffix("\","))
        .find(|line| line.split_whitespace().next() == Some(dependency))?;
    if let Some(version) = reference.split_whitespace().nth(1) {
        return Some(version.into());
    }
    let name = format!("name = \"{dependency}\"");
    let mut matches = blocks
        .iter()
        .filter(|block| block.lines().any(|line| line.trim() == name));
    let block = matches.next()?;
    if matches.next().is_some() {
        return None;
    }
    block.lines().find_map(|line| {
        line.trim()
            .strip_prefix("version = \"")?
            .strip_suffix('"')
            .map(str::to_owned)
    })
}

fn raster_metadata(
    source_hash: &str,
    source_file: Option<&str>,
    output: &str,
    format: &str,
    width: u32,
    height: u32,
    font_faces: usize,
    used_fonts: &[SvgUsedFontIdentity],
) -> Result<Value, EngineError> {
    let metadata = json!({
        "format": format, "output_path": output, "source_svg_path": source_file,
        "source_svg_sha256": source_hash,
        "backend": "in-process resvg via audited crate::svg_png",
        "resvg_version": root_dependency_version("resvg"),
        "backend_lockfile_sha256": sha256_hex_bytes(LOCKFILE.as_bytes()),
        "scale": 1.0, "drop_dotplot_metadata": false, "width": width, "height": height,
        "font_face_count": font_faces,
        "used_font_face_count": used_fonts.len(),
        "font_identities": used_fonts,
        "font_identity_status": FONT_IDENTITY_STATUS,
        "font_digest_convention": FONT_DIGEST_CONVENTION,
        "font_audit_scope": "selected positioned-glyph layout fonts, including fallback and nested SVG trees; not proof of visible pixels or complete glyph coverage",
        "font_loading_policy": "system fonts plus GENTLE_SVG_FONT_FILE and GENTLE_SVG_FONT_DIR; generic family overrides honored by svg_png",
        "font_reproducibility": "recorded_not_enforced_or_replayed",
    });
    validate_raster_font_metadata(&metadata)?;
    Ok(metadata)
}

fn validate_raster_font_metadata(metadata: &Value) -> Result<(), EngineError> {
    let identities = metadata
        .get("font_identities")
        .and_then(Value::as_array)
        .ok_or_else(|| invalid("raster metadata requires actual glyph-used font identities"))?;
    let used_count = identities.len() as u64;
    if identities.is_empty()
        || identities.len() > 512
        || metadata.get("used_font_face_count").and_then(Value::as_u64) != Some(used_count)
        || metadata
            .get("font_face_count")
            .and_then(Value::as_u64)
            .is_none_or(|count| count < used_count)
        || metadata.get("font_identity_status").and_then(Value::as_str)
            != Some(FONT_IDENTITY_STATUS)
        || metadata
            .get("font_digest_convention")
            .and_then(Value::as_str)
            != Some(FONT_DIGEST_CONVENTION)
        || metadata.get("font_reproducibility").and_then(Value::as_str)
            != Some("recorded_not_enforced_or_replayed")
    {
        return Err(invalid(
            "missing or inconsistent glyph-used font audit for a text-bearing TSS raster",
        ));
    }
    let mut seen = BTreeSet::new();
    for identity in identities {
        let fields = identity
            .as_object()
            .ok_or_else(|| invalid("invalid used-font identity"))?;
        if fields.len() != 4
            || !fields.keys().all(|key| {
                matches!(
                    key.as_str(),
                    "families" | "post_script_name" | "face_index" | "sha256"
                )
            })
            || identity
                .get("post_script_name")
                .and_then(Value::as_str)
                .is_none()
        {
            return Err(invalid(
                "used-font identities require families, post_script_name, face_index and sha256 only",
            ));
        }
        let families = identity
            .get("families")
            .and_then(Value::as_array)
            .ok_or_else(|| invalid("used-font families must be an array"))?;
        if families.iter().any(|family| !family.is_string()) {
            return Err(invalid("used-font family names must be strings"));
        }
        let hash = identity
            .get("sha256")
            .and_then(Value::as_str)
            .ok_or_else(|| invalid("used-font identity lacks a source digest"))?;
        digest_field(hash, "used-font source digest")?;
        let index = identity
            .get("face_index")
            .and_then(Value::as_u64)
            .filter(|index| u32::try_from(*index).is_ok())
            .ok_or_else(|| invalid("used-font face index must be a nonnegative u32"))?;
        if !seen.insert((hash, index)) {
            return Err(invalid("duplicate used-font source/face identity"));
        }
    }
    Ok(())
}

fn write_page(
    inventory: &mut Inventory<'_>,
    stem: &str,
    page: &TssRenderedPage,
    formats: &[TssExportFormat],
    metadata: &mut Vec<Value>,
    should_continue: &mut dyn FnMut() -> bool,
) -> Result<PageIndex, EngineError> {
    let prefix = format!("{stem}.page-{:04}", page.page_number);
    let svg_name = format!("{prefix}.svg");
    let svg_hash = sha256_hex_bytes(page.svg.as_bytes());
    let keep_svg = formats.contains(&TssExportFormat::Svg);
    let source_file = keep_svg.then_some(svg_name.as_str());
    let mut files = Vec::new();
    if keep_svg {
        checkpoint(should_continue)?;
        inventory.bytes(&svg_name, page.svg.as_bytes())?;
        let (width, height) = svg_dimensions(&page.svg)?;
        metadata.push(json!({
            "format": "svg", "output_path": svg_name, "width": width, "height": height,
            "backend": "gentle_render::tss_profiles::render_tss_profile_pages",
            "exporter_revision": option_env!("GENTLE_SOURCE_REVISION").unwrap_or("unknown"),
            "font_resolution": "deferred_to_svg_viewer", "font_identities": null,
            "font_identity_status": "SVG retains renderer font-family declarations; viewer-selected fonts are not pinned",
        }));
        files.push(svg_name.clone());
    }
    let options = SvgPngRenderOptions::default();
    if formats.contains(&TssExportFormat::Png) {
        checkpoint(should_continue)?;
        let (png, used_fonts) = crate::svg_png::render_svg_to_png_bytes_audited(&page.svg, options)
            .map_err(|e| invalid(format!("PNG rendering failed: {e}")))?;
        let name = format!("{prefix}.png");
        let entry = raster_metadata(
            &svg_hash,
            source_file,
            &name,
            "png",
            png.width,
            png.height,
            png.font_face_count,
            &used_fonts,
        )?;
        inventory.bytes(&name, &png.bytes)?;
        metadata.push(entry);
        files.push(name);
    }
    if formats.contains(&TssExportFormat::Pdf) {
        checkpoint(should_continue)?;
        // The existing PDF helper takes a filename. The temporary SVG is never
        // published or included as an alleged output when SVG was not requested.
        let mut source = tempfile::Builder::new()
            .prefix(".render-")
            .suffix(".svg")
            .tempfile_in(inventory.directory)
            .map_err(|e| io_error("create temporary PDF source", e))?;
        source
            .write_all(page.svg.as_bytes())
            .map_err(|e| io_error("write temporary PDF source", e))?;
        source
            .flush()
            .map_err(|e| io_error("flush temporary PDF source", e))?;
        let name = format!("{prefix}.pdf");
        let (path, file) = inventory.reserve(&name)?;
        drop(file);
        let (summary, used_fonts) =
            crate::svg_pdf::render_svg_file_to_pdf_audited(source.path(), &path, options)
                .map_err(|e| invalid(format!("raster-backed PDF rendering failed: {e}")))?;
        File::open(&path)
            .and_then(|file| file.sync_all())
            .map_err(|e| io_error("sync PDF", e))?;
        inventory.record(&name, &path)?;
        let mut entry = raster_metadata(
            &svg_hash,
            source_file,
            &name,
            "pdf",
            summary.width,
            summary.height,
            summary.font_face_count,
            &used_fonts,
        )?;
        entry["pdf_representation"] = json!("single-page raster-backed RGB image; not vector PDF");
        entry["pdf_image_encoding"] = json!("FlateDecode (lossless zlib-compressed RGB)");
        entry["pdf_helper"] = json!("crate::svg_pdf::render_svg_file_to_pdf_audited");
        entry["page_width_pt"] = json!(summary.page_width_pt);
        entry["page_height_pt"] = json!(summary.page_height_pt);
        entry["uri_link_count"] = json!(summary.uri_link_count);
        metadata.push(entry);
        files.push(name);
    }
    Ok(PageIndex {
        page_number: page.page_number,
        page_count: page.page_count,
        promoter_ids: page.promoter_ids.clone(),
        files,
    })
}

fn readme(report: &TssProfileReport) -> Result<String, EngineError> {
    let source = report.source.as_ref();
    let assessed = |value: Option<&str>| value.map(tsv).unwrap_or_else(|| "not_assessed".into());
    let provenance = format!(
        "## Source Provenance\n\n\
Input manifest SHA-256: {}\n\n\
Input source schema: {}\n\n\
Input source revision: {}\n\n\
Input dataset ID: {}\n\n\
Bundle-source producer SHA-256 (source-declared): {}\n\n\
Separately declared annotation release: {}\n\n\
Profile producer executable SHA-256: {}\n\n\
The named receipt input_manifest_sha256 is copied from the report's explicit bundle_manifest\n\
input binding (legacy manifest role also accepted), which binds the original manifest bytes.\n\
When report.source is present its manifest_sha256 must match that binding exactly; its\n\
source_revision is retained separately from the scoring producer and exporter build revisions.\n\
A legacy report without source metadata still binds its input manifest, but source_revision\n\
remains null. Missing optional provenance is not_assessed (JSON/TSV null), not a guessed value.\n\
In particular, annotation release is never extracted from words in a genome label.\n\
The optional producer_executable receipt input binds the supplied profile producer binary\n\
digest; bundle_source_producer binds the separate source-declared producer digest when present.\n\
Neither substitutes for receipt.executable_sha256, which hashes the running exporter. These\n\
identities are copied from the report, not independently reverified against original input files.\n\n",
        manifest_digest(&report.inputs)?,
        assessed(source.map(|value| value.schema.as_str())),
        assessed(source_revision(report)),
        assessed(source.and_then(|value| value.dataset_id.as_deref())),
        assessed(source.and_then(|value| value.producer_sha256.as_deref())),
        assessed(report.reference.annotation_release.as_deref()),
        assessed(report.producer_executable_sha256.as_deref()),
    );
    let mut descriptions = String::from(
        "\n## Selected TSS Descriptions\n\nSelection labels, legends, criteria and optional factors below are descriptive report metadata,\nnot an independent evidence assessment. Selection does not establish TSS usage, direct binding\nor promoter activity. A selected TSS without recorded evidence gets only a generic label and\nlegend; its criterion and factor remain not_assessed. No CUT&RUN support is inferred.\n\n",
    );
    for window in report.windows.iter().filter(|window| window.selected) {
        let [label, legend, criterion, factor] = selection_fields(window);
        descriptions.push_str(&format!(
            "### {} / {} / {}\n\nSelected-panel label: {}\n\nLegend: {}\n\nRecorded criterion: {}\n\nRecorded factor: {}\n\n",
            tsv(&window.record.gene_symbol), tsv(&window.record.gene_id), tsv(&window.record.promoter_id),
            assessed(label), assessed(legend), assessed(criterion), assessed(factor),
        ));
    }
    if !report.windows.iter().any(|window| window.selected) {
        descriptions.push_str("No TSS is marked selected in the supplied report.\n");
    }
    Ok(format!(
        "# TSS TFBS Profile Export\n\n\
This export projects the supplied report without rescoring or reading source sequences.\n\n\
## Interpretation Limits\n\n{NON_CLAIMS}\n\n\
Bundle consistency and hashes establish internal consistency, not independent reference authenticity.\n\
The report's verification status is retained verbatim in report.json. Nothing here adds biological validation.\n\n\
{provenance}\
## Files And Coordinates\n\n\
report.json is the complete computational report. Gene JSON files are protocol-compatible subsets,\n\
including their comparisons, original raw scores, exact panel, normalization metadata, source\n\
provenance, optional producer binary digest, selection evidence and input bindings.\n\
Each resolved_matrix digest is lowercase SHA-256 of serde_json::to_vec of the tuple\n\
(accession, Some(exact factor name), matrix_counts), serialized as [accession,name,counts].\n\
Counts retain A/C/G/T column order and f64 JSON serialization; this is not a counts-only digest.\n\
index.json lists genes, TSS memberships, pages and data-file SHA-256 values. Gene order and TSS order\n\
follow the supplied report; filenames combine safe symbol/ID fragments and a full identity digest.\n\
The index repeats named source provenance and selected-panel descriptions for cross-gene inspection.\n\
comparisons.tsv is the cross-gene table; each gene also has its own comparison TSV.\n\n\
Score TSV is long-form: one row per sequence-base window start, matrix and local motif strand.\n\
local_window_start_0based is a start, never a motif center or an inferred motif 5-prime endpoint.\n\
tss_relative_window_start_bp and genomic_window_start_1based use protocol TssGeometry methods.\n\
The common start coordinate is identical for both local motif strands. Genomic positions decrease\n\
along a minus-transcript window. Motif intervals are inclusive, 1-based and stored ascending;\n\
motif_5prime_genomic_1based is the opposite interval endpoint for the local reverse motif.\n\
Both local and genomic motif strands are explicit. An unavailable raw score is literal null, not zero.\n\
Non-fitting starts (including the trailing motif-length-minus-one bases) have null scores and motif\n\
interval/endpoint fields, with availability=outside_window. Internal nulls remain unavailable.\n\
display_score alone applies panel.clip_negative (max(raw_score, 0)); raw_score is unchanged.\n\
Score units and all computation policies are in the panel and score_policy, not inferred from heights.\n\
Comparison values are copied, not recomputed: their method retains strand, smoothing and raw-input\n\
policy; null correlations carry undefined_reason. Display clipping does not alter these values.\n\
Every unordered same-factor matrix pair must have one comparison on each local strand. Paired\n\
counts must match starts where both raw vectors are available; excluded counts cover the longer\n\
vector's start domain minus paired starts. Validation checks these masks, not correlation formulas.\n\
Every TSV begins with # comment lines carrying canonical non-claims and coordinate/method policies.\n\
Provenance and selected-panel legends are also included in those comments. Selection columns\n\
preserve the descriptive label, legend, criterion and factor; unselected rows use null for all four.\n\
Skip these lines before reading the column header. TSV text escapes backslash, tab, CR and LF\n\
as \\\\, \\t, \\r and \\n respectively; a leading # in a text cell is escaped as \\#.\n\n\
## Rendering And Replay\n\n\
SVG pages use gentle_render::tss_profiles::render_tss_profile_pages. Optional PNG uses in-process\n\
resvg at scale 1 with metadata stripping disabled. PDF is a single-page raster-backed RGB image,\n\
not a vector PDF. Its FlateDecode/zlib compression preserves every RGB pixel and the resolution.\n\
Receipt metadata records that encoding, actual dimensions, available font-face counts and the\n\
root's locked resvg version. Audited PNG/PDF helpers inspect positioned glyphs in the same parsed\n\
rendering tree, including fallback and nested SVG trees. Each font_identities entry records its\n\
family names, PostScript name, face_index and lowercase SHA-256 of the complete font source/container\n\
bytes. The face index is recorded separately, not mixed into that digest. used_font_face_count counts\n\
these distinct used source/face identities; font_face_count counts available faces, not used faces.\n\
An unbound glyph-used face fails export instead of leaving a partial audit. This is selected\n\
layout-font evidence, not proof that every glyph produced visible pixels or that glyph coverage is\n\
complete. Font files and host paths are not copied into the export. For font-matched replay, supply\n\
the exact recorded font bytes and compare the newly selected identities, backend and options.\n\
Font pinning is not enforced; cross-host byte-identical replay is NOT verified by recording these\n\
identities alone. SVG-only pages retain font-family declarations and viewer-dependent font selection;\n\
their font_identities remain null. PNG/PDF font bindings do not pin how another viewer renders SVG.\n\n\
For report-only replay, load report.json as TssProfileReport and export-request.json as\n\
ExportTssProfilesRequest. Replace output_dir (stored portably as '.') with a fresh output directory,\n\
then invoke export_tss_profiles through the shared export operation. No scores are recomputed.\n\
To recompute instead, separately recover the exact hash-bound manifest, FASTA, panel, optional\n\
selection and registry sources. The original inputs are not copied into this export. Use the\n\
producer revision and producer lockfile digest in report.json, not the exporter revision as a substitute.\n\n\
Producer revision: {}\n\
Exporter revision: {}\n\
Receipt executable_sha256 hashes the running exporter executable. Receipt lockfile_sha256 hashes\n\
the exporter build's embedded Cargo.lock; index.json separately retains the producer lockfile digest.\n\
No timestamps or staging/host paths enter the result identity. report_sha256 hashes the exact\n\
compact report.json bytes. The receipt binds all published files except receipt.json itself.\n\
The index hashes neither itself nor the receipt, avoiding a circular hash. Verification checks\n\
the complete inventory, file hashes and report bindings; it is not a digital signature or proof\n\
of authenticity against an attacker who can rewrite both files and the receipt.\n\
{descriptions}",
        tsv(&report.producer_revision),
        option_env!("GENTLE_SOURCE_REVISION").unwrap_or("unknown")
    ))
}

/// Export an already computed report into a fresh directory, without rescoring.
pub fn export_tss_profiles(
    report: &TssProfileReport,
    request: &ExportTssProfilesRequest,
) -> Result<TssProfileReceipt, EngineError> {
    export_tss_profiles_with_cancel(report, request, &mut || true)
}

/// Cancellable export with checks between report/render/file stages and TSV chunks.
///
/// A false callback abandons the staging directory. The existing renderer and
/// raster helpers are synchronous, so a single page conversion is not interruptible.
/// The final publish is a sibling directory rename after a second freshness check;
/// this is not a defense against a hostile process replacing parent directories.
pub fn export_tss_profiles_with_cancel(
    report: &TssProfileReport,
    request: &ExportTssProfilesRequest,
    should_continue: &mut dyn FnMut() -> bool,
) -> Result<TssProfileReceipt, EngineError> {
    checkpoint(should_continue)?;
    if request.context_manifest.is_some() {
        return Err(invalid(
            "Resolve context_manifest through the shared TSS export operation before rendering",
        ));
    }
    preflight_tss_export(request)?;
    validate_tss_profile_report(report)?;
    validate_rendering(report, request)?;
    let genes = genes(report)?;
    let inputs = input_bindings(report)?;
    checkpoint(should_continue)?;
    let pages = render_tss_profile_pages(report, &request.rendering)
        .map_err(|e| invalid(format!("page rendering failed: {e}")))?;
    validate_pages(
        &pages,
        &genes,
        request.formats.iter().any(|f| *f != TssExportFormat::Svg),
    )?;
    checkpoint(should_continue)?;
    let executable =
        std::env::current_exe().map_err(|e| io_error("locate exporter executable", e))?;
    let executable_sha256 =
        sha256_file_hex(&executable).map_err(|e| io_error("hash exporter executable", e))?;
    let output = destination(&request.output_dir)?;
    let parent = output
        .parent()
        .ok_or_else(|| invalid("output has no parent directory"))?;
    let staging = tempfile::Builder::new()
        .prefix(".gentle-tss-export-")
        .tempdir_in(parent)
        .map_err(|e| io_error("create fresh staging directory", e))?;
    let mut inventory = Inventory {
        directory: staging.path(),
        hashes: BTreeMap::new(),
        bytes: 0,
    };
    inventory.json(REPORT_FILE, report)?;
    let report_sha256 = inventory.hashes[REPORT_FILE].clone();
    let portable_request = ExportTssProfilesRequest {
        context_manifest: None,
        output_dir: ".".into(),
        rendering: request.rendering.clone(),
        formats: request.formats.clone(),
    };
    inventory.json(REQUEST_FILE, &portable_request)?;
    inventory.bytes(README_FILE, readme(report)?.as_bytes())?;
    inventory.write(COMPARISON_FILE, |writer| {
        write_comparisons(writer, report, &report.windows)
    })?;
    let mut gene_index = Vec::new();
    let mut render_metadata = Vec::new();
    for gene in &genes {
        checkpoint(should_continue)?;
        let report_file = format!("{}.json", gene.stem);
        let scores_file = format!("{}.scores.tsv", gene.stem);
        let comparisons_file = format!("{}.comparisons.tsv", gene.stem);
        inventory.json(&report_file, &GeneReport::new(report, gene))?;
        inventory.write(&scores_file, |writer| {
            write_scores(writer, report, &gene.windows, should_continue)
        })?;
        inventory.write(&comparisons_file, |writer| {
            write_comparisons(writer, report, gene.windows.iter().copied())
        })?;
        let mut gene_pages = Vec::new();
        for page in pages.iter().filter(|page| page.gene_id == gene.id) {
            gene_pages.push(write_page(
                &mut inventory,
                &gene.stem,
                page,
                &request.formats,
                &mut render_metadata,
                should_continue,
            )?);
        }
        gene_index.push(GeneIndex {
            gene_id: gene.id.into(),
            gene_symbol: gene.symbol.into(),
            promoter_ids: gene
                .windows
                .iter()
                .map(|w| w.record.promoter_id.clone())
                .collect(),
            selections: selection_index(gene),
            report: report_file,
            scores: scores_file,
            comparisons: comparisons_file,
            pages: gene_pages,
        });
    }
    checkpoint(should_continue)?;
    let index = ExportIndex {
        schema: INDEX_SCHEMA.into(),
        report_sha256: report_sha256.clone(),
        input_manifest_sha256: manifest_digest(&report.inputs)?.into(),
        source_revision: source_revision(report).map(str::to_owned),
        source: report.source.clone(),
        reference: report.reference.clone(),
        producer_revision: report.producer_revision.clone(),
        producer_executable_sha256: report.producer_executable_sha256.clone(),
        producer_lockfile_sha256: report.lockfile_sha256.clone(),
        score_policy: report.score_policy.clone(),
        rendering: request.rendering.clone(),
        formats: request.formats.clone(),
        tss_count: report.windows.len(),
        page_count: pages.len(),
        genes: gene_index,
        outputs: inventory.hashes.clone(),
    };
    bounded_json(&index, MAX_METADATA_BYTES)?;
    inventory.json(INDEX_FILE, &index)?;
    let receipt = TssProfileReceipt {
        schema: RECEIPT_SCHEMA.into(),
        report_sha256,
        input_manifest_sha256: index.input_manifest_sha256.clone(),
        source_revision: index.source_revision.clone(),
        producer_revision: report.producer_revision.clone(),
        exporter_revision: option_env!("GENTLE_SOURCE_REVISION")
            .unwrap_or("unknown")
            .into(),
        executable_sha256,
        lockfile_sha256: sha256_hex_bytes(LOCKFILE.as_bytes()),
        inputs,
        rendering: request.rendering.clone(),
        renderer: "gentle_render::tss_profiles::render_tss_profile_pages".into(),
        render_metadata,
        outputs: inventory.hashes.clone(),
        tss_count: report.windows.len(),
        page_count: pages.len(),
        non_claims: NON_CLAIMS.into(),
    };
    bounded_json(&receipt, MAX_METADATA_BYTES)?;
    inventory.json(RECEIPT_FILE, &receipt)?;
    checkpoint(should_continue)?;
    verify_tss_profile_receipt(staging.path(), &receipt)?;
    checkpoint(should_continue)?;
    drop(inventory);
    checked_directory(parent)?;
    check_fresh_destination(&output)?;
    // Remove only a verified empty directory, using rmdir (never recursive
    // removal). rename also refuses a destination that became nonempty.
    if output.exists() {
        fs::remove_dir(&output)
            .map_err(|e| io_error("remove empty destination before publish", e))?;
    }
    fs::rename(staging.path(), &output).map_err(|e| io_error("publish completed export", e))?;
    let _ = staging.keep();
    Ok(receipt)
}

/// Validate receipt structure without touching the filesystem. Empty inventories
/// and self-hashes are invalid; a receipt must bind a complete nonempty export.
pub fn validate_tss_profile_receipt(receipt: &TssProfileReceipt) -> Result<(), EngineError> {
    if receipt.schema != RECEIPT_SCHEMA
        || receipt.non_claims != NON_CLAIMS
        || receipt.outputs.len() > MAX_OUTPUT_FILES
        || receipt.outputs.contains_key(RECEIPT_FILE)
        || receipt.tss_count == 0
        || receipt.tss_count > MAX_WINDOWS
        || receipt.page_count == 0
        || receipt.page_count > MAX_OUTPUT_FILES / 3
        || receipt.inputs.is_empty()
        || receipt.render_metadata.is_empty()
        || receipt.render_metadata.len() > MAX_OUTPUT_FILES
    {
        return Err(invalid(
            "invalid or empty receipt, self-hash, counts or non-claims",
        ));
    }
    for name in [
        REPORT_FILE,
        INDEX_FILE,
        REQUEST_FILE,
        README_FILE,
        COMPARISON_FILE,
    ] {
        if !receipt.outputs.contains_key(name) {
            return Err(invalid(format!("receipt omits required output {name}")));
        }
    }
    for hash in [
        &receipt.report_sha256,
        &receipt.input_manifest_sha256,
        &receipt.executable_sha256,
        &receipt.lockfile_sha256,
    ] {
        digest_field(hash, "receipt digest")?;
    }
    if receipt.outputs.get(REPORT_FILE) != Some(&receipt.report_sha256) {
        return Err(invalid(
            "receipt report digest disagrees with its output inventory",
        ));
    }
    for (name, hash) in &receipt.outputs {
        if !portable_name(name) {
            return Err(invalid("receipt contains an unsafe output filename"));
        }
        digest_field(hash, "output digest")?;
    }
    for revision in [&receipt.producer_revision, &receipt.exporter_revision] {
        text_field(revision, "build revision")?;
    }
    if let Some(revision) = &receipt.source_revision {
        provenance_text(revision, "receipt source revision")?;
    }
    if receipt.renderer != "gentle_render::tss_profiles::render_tss_profile_pages" {
        return Err(invalid("unknown receipt renderer"));
    }
    validate_bindings(&receipt.inputs)?;
    if manifest_digest(&receipt.inputs)? != receipt.input_manifest_sha256 {
        return Err(invalid(
            "receipt input_manifest_sha256 disagrees with its input bindings",
        ));
    }
    for value in &receipt.render_metadata {
        validate_metadata(value)?;
        if matches!(
            value.get("format").and_then(Value::as_str),
            Some("png" | "pdf")
        ) {
            validate_raster_font_metadata(value)?;
        }
    }
    bounded_json(receipt, MAX_METADATA_BYTES)
}

fn read_json<T: serde::de::DeserializeOwned>(path: &Path, limit: u64) -> Result<T, EngineError> {
    let metadata = fs::symlink_metadata(path).map_err(|e| io_error("inspect audit JSON", e))?;
    if metadata.file_type().is_symlink() || !metadata.is_file() || metadata.len() > limit {
        return Err(invalid(
            "audit JSON must be a bounded ordinary file, not a symlink",
        ));
    }
    let file = File::open(path).map_err(|e| io_error("open audit JSON", e))?;
    serde_json::from_reader(BufReader::new(file))
        .map_err(|e| invalid(format!("invalid audit JSON: {e}")))
}

fn same_json(left: &impl Serialize, right: &impl Serialize) -> Result<bool, EngineError> {
    Ok(
        serde_json::to_vec(left).map_err(|e| invalid(e.to_string()))?
            == serde_json::to_vec(right).map_err(|e| invalid(e.to_string()))?,
    )
}

/// Verify actual file hashes, the exact inventory and report/index/receipt links.
///
/// The caller may supply a separately trusted receipt to detect changes to its
/// on-disk copy as well. Without an externally trusted digest or signature this
/// proves internal consistency, not producer authenticity. No files are modified.
pub fn verify_tss_profile_receipt(
    output_dir: &Path,
    receipt: &TssProfileReceipt,
) -> Result<(), EngineError> {
    validate_tss_profile_receipt(receipt)?;
    let output = checked_directory(output_dir)?;
    let on_disk: TssProfileReceipt = read_json(&output.join(RECEIPT_FILE), MAX_METADATA_BYTES)?;
    if !same_json(receipt, &on_disk)? {
        return Err(invalid("receipt.json differs from the supplied receipt"));
    }
    let mut actual = BTreeSet::new();
    let mut total = 0u64;
    for entry in fs::read_dir(&output).map_err(|e| io_error("list export inventory", e))? {
        let entry = entry.map_err(|e| io_error("read export entry", e))?;
        let name = entry
            .file_name()
            .into_string()
            .map_err(|_| invalid("non-UTF-8 export filename"))?;
        if !portable_name(&name) || (name != RECEIPT_FILE && !receipt.outputs.contains_key(&name)) {
            return Err(invalid("export contains an unexpected or unsafe artifact"));
        }
        let metadata =
            fs::symlink_metadata(entry.path()).map_err(|e| io_error("inspect export entry", e))?;
        if metadata.file_type().is_symlink() || !metadata.is_file() {
            return Err(invalid(
                "export inventory contains a symlink, subdirectory or special file",
            ));
        }
        total = total
            .checked_add(metadata.len())
            .filter(|n| *n <= MAX_EXPORT_BYTES)
            .ok_or_else(|| invalid("export inventory exceeds its byte budget"))?;
        actual.insert(name);
    }
    if actual.len() != receipt.outputs.len() + 1 {
        return Err(invalid("export is missing an inventoried artifact"));
    }
    for (name, expected) in &receipt.outputs {
        let actual = sha256_file_hex(&output.join(name))
            .map_err(|e| io_error("hash exported artifact", e))?;
        if actual != *expected {
            return Err(invalid(format!("SHA-256 mismatch for {name}")));
        }
    }
    let report: TssProfileReport = read_json(&output.join(REPORT_FILE), MAX_REPORT_BYTES)?;
    validate_tss_profile_report(&report)?;
    let request: ExportTssProfilesRequest =
        read_json(&output.join(REQUEST_FILE), MAX_METADATA_BYTES)?;
    validate_rendering(&report, &request)?;
    if request.output_dir != "."
        || !same_json(&request.rendering, &receipt.rendering)?
        || receipt.input_manifest_sha256 != manifest_digest(&report.inputs)?
        || receipt.source_revision.as_deref() != source_revision(&report)
        || report.producer_revision != receipt.producer_revision
        || report.windows.len() != receipt.tss_count
        || !same_json(&input_bindings(&report)?, &receipt.inputs)?
    {
        return Err(invalid(
            "receipt inputs, source provenance, counts, producer or rendering disagree with the source report/request",
        ));
    }
    let index: ExportIndex = read_json(&output.join(INDEX_FILE), MAX_METADATA_BYTES)?;
    let mut indexed_hashes = receipt.outputs.clone();
    indexed_hashes.remove(INDEX_FILE);
    if index.schema != INDEX_SCHEMA
        || index.report_sha256 != receipt.report_sha256
        || index.input_manifest_sha256 != receipt.input_manifest_sha256
        || index.source_revision != receipt.source_revision
        || !same_json(&index.source, &report.source)?
        || index.producer_revision != report.producer_revision
        || index.producer_executable_sha256 != report.producer_executable_sha256
        || index.producer_lockfile_sha256 != report.lockfile_sha256
        || index.reference != report.reference
        || index.score_policy != report.score_policy
        || index.tss_count != receipt.tss_count
        || index.page_count != receipt.page_count
        || !same_json(&index.rendering, &receipt.rendering)?
        || index.formats != request.formats
        || index.outputs != indexed_hashes
    {
        return Err(invalid(
            "index disagrees with the bound report, policies or data inventory",
        ));
    }
    let genes = genes(&report)?;
    if index.genes.len() != genes.len() {
        return Err(invalid("index gene count mismatch"));
    }
    let mut expected = [
        REPORT_FILE,
        INDEX_FILE,
        REQUEST_FILE,
        README_FILE,
        COMPARISON_FILE,
    ]
    .into_iter()
    .map(str::to_owned)
    .collect::<BTreeSet<_>>();
    let mut image_files = BTreeSet::new();
    let mut page_count = 0usize;
    for (entry, gene) in index.genes.iter().zip(&genes) {
        let promoters: Vec<_> = gene
            .windows
            .iter()
            .map(|w| w.record.promoter_id.as_str())
            .collect();
        if entry.gene_id != gene.id
            || entry.gene_symbol != gene.symbol
            || entry
                .promoter_ids
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>()
                != promoters
            || !same_json(&entry.selections, &selection_index(gene))?
            || entry.report != format!("{}.json", gene.stem)
            || entry.scores != format!("{}.scores.tsv", gene.stem)
            || entry.comparisons != format!("{}.comparisons.tsv", gene.stem)
            || entry.pages.is_empty()
        {
            return Err(invalid(
                "index gene order, membership, selection evidence or filenames disagree with the report",
            ));
        }
        expected.extend([
            entry.report.clone(),
            entry.scores.clone(),
            entry.comparisons.clone(),
        ]);
        let subset: TssProfileReport = read_json(&output.join(&entry.report), MAX_REPORT_BYTES)?;
        if !same_json(&subset, &GeneReport::new(&report, gene))? {
            return Err(invalid(
                "per-gene JSON is not an unchanged subset of the source report",
            ));
        }
        let mut seen_promoters = BTreeSet::new();
        for (page_index, page) in entry.pages.iter().enumerate() {
            page_count += 1;
            if page.page_number != page_index + 1
                || page.page_count != entry.pages.len()
                || page.promoter_ids.is_empty()
            {
                return Err(invalid("invalid index pagination"));
            }
            let mut members = BTreeSet::new();
            for promoter in &page.promoter_ids {
                if !members.insert(promoter) || !promoters.contains(&promoter.as_str()) {
                    return Err(invalid("invalid index page membership"));
                }
                seen_promoters.insert(promoter.as_str());
            }
            let files = [
                (TssExportFormat::Svg, "svg"),
                (TssExportFormat::Png, "png"),
                (TssExportFormat::Pdf, "pdf"),
            ]
            .into_iter()
            .filter(|(format, _)| request.formats.contains(format))
            .map(|(_, extension)| format!("{}.page-{:04}.{extension}", gene.stem, page.page_number))
            .collect::<Vec<_>>();
            if page.files != files {
                return Err(invalid(
                    "index page formats/filenames disagree with request",
                ));
            }
            expected.extend(files.iter().cloned());
            image_files.extend(files);
        }
        if promoters.iter().any(|id| !seen_promoters.contains(id)) {
            return Err(invalid("index pages omit a TSS"));
        }
    }
    if page_count != receipt.page_count || expected != receipt.outputs.keys().cloned().collect() {
        return Err(invalid(
            "receipt output/page inventory is incomplete or contains extra files",
        ));
    }
    let mut described_images = BTreeSet::new();
    for metadata in &receipt.render_metadata {
        let name = metadata
            .get("output_path")
            .and_then(Value::as_str)
            .ok_or_else(|| invalid("render metadata lacks output_path"))?;
        if !image_files.contains(name) || !described_images.insert(name.to_owned()) {
            return Err(invalid(
                "render metadata contains duplicate or unknown outputs",
            ));
        }
        if metadata.get("format").and_then(Value::as_str) != name.rsplit('.').next() {
            return Err(invalid("render metadata format does not match its file"));
        }
    }
    if described_images != image_files {
        return Err(invalid("render metadata omits output images"));
    }
    Ok(())
}

/// Load the on-disk receipt and verify its complete export without modifying it.
pub fn read_and_verify_tss_profile_receipt(
    output_dir: &Path,
) -> Result<TssProfileReceipt, EngineError> {
    let output = checked_directory(output_dir)?;
    let receipt = read_json(&output.join(RECEIPT_FILE), MAX_METADATA_BYTES)?;
    verify_tss_profile_receipt(&output, &receipt)?;
    Ok(receipt)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    // Entirely hand-crafted export-contract fixture, recreated by this function.
    // These invented matrices, scores and reference labels are NOT JASPAR data
    // or scorer acceptance evidence. They exercise serialization, coordinates,
    // rendering and receipts only; no external test_files/tests fixtures are used.
    pub(crate) fn synthetic_report() -> TssProfileReport {
        let factors = [
            ("MA0001.1", "SYNTH_A", "Synthetic A version 1"),
            ("MA0001.2", "SYNTH_A", "Synthetic A version 2"),
            ("MA0002.1", "SYNTH_B", "Synthetic B"),
        ]
        .into_iter()
        .enumerate()
        .map(|(i, (accession, factor, label))| JasparPanelTrack {
            track_id: None,
            provider_kind: None,
            factor_label: None,
            source_id: accession.into(),
            factor_id: factor.into(),
            label: label.into(),
            display_order: i + 1,
            color_hint: None,
            score_kind: None,
        })
        .collect::<Vec<_>>();
        let panel = JasparTargetPanel {
            schema: PANEL_SCHEMA.into(),
            panel_id: "synthetic-export-only".into(),
            label: "Synthetic export contract fixture".into(),
            score_kind: "llr_bits".into(),
            clip_negative: true,
            scale_mode: TssScaleMode::Independent,
            strand_policy: TssStrandPolicy::Both,
            calibration_state: TssCalibrationState::MatrixSpecific,
            calibration_statement: "Synthetic model outputs, not calibrated between matrices"
                .into(),
            calibration_id: None,
            calibration_sha256: None,
            top_hit_count: 2,
            factors,
        };
        let matrices = panel
            .factors
            .iter()
            .enumerate()
            .map(|(i, specification)| {
                let counts = vec![[8.0, 1.0, 1.0, 2.0]; 3 - i];
                ResolvedTssMatrix {
                    specification: specification.clone(),
                    version: specification.source_id.rsplit('.').next().unwrap().into(),
                    consensus: "A".repeat(counts.len()),
                    matrix_sha256: sha256_hex_bytes(
                        &serde_json::to_vec(&(
                            &specification.source_id,
                            &Some(specification.factor_id.clone()),
                            &counts,
                        ))
                        .unwrap(),
                    ),
                    matrix_counts: counts,
                }
            })
            .collect::<Vec<_>>();
        let binding = |role: &str, name: &str, bytes: &[u8]| TssInputBinding {
            role: role.into(),
            name: name.into(),
            sha256: sha256_hex_bytes(bytes),
        };
        let panel_sha256 = sha256_hex_bytes(&serde_json::to_vec(&panel).unwrap());
        let mut inputs = vec![
            binding(
                "manifest",
                "manifest.json",
                b"hand-crafted synthetic manifest",
            ),
            binding("fasta", "windows.fa", b"hand-crafted synthetic FASTA"),
        ];
        inputs.push(TssInputBinding {
            role: "panel".into(),
            name: "panel.json".into(),
            sha256: panel_sha256.clone(),
        });
        let signals = [
            (
                vec![Some(-2.0), None, Some(3.0)],
                vec![Some(4.0), None, Some(-1.0)],
            ),
            (
                vec![Some(1.0), None, Some(2.0), Some(5.0)],
                vec![Some(-3.0), None, Some(1.0), Some(2.0)],
            ),
            (
                vec![Some(-4.0), Some(-3.0), Some(-2.0), Some(-1.0), Some(0.0)],
                vec![Some(5.0), Some(4.0), Some(3.0), Some(2.0), Some(1.0)],
            ),
        ];
        let peaks = |scores: &[Option<f64>]| {
            let mut result = scores
                .iter()
                .enumerate()
                .filter_map(|(i, value)| {
                    value.map(|score| TssPeak {
                        local_start_0based: i,
                        score,
                    })
                })
                .collect::<Vec<_>>();
            result.sort_by(|a, b| {
                b.score
                    .total_cmp(&a.score)
                    .then(a.local_start_0based.cmp(&b.local_start_0based))
            });
            result.truncate(2);
            result
        };
        let windows = [(TssStrand::Plus, "gene-a", "SYN/A", "promoter-plus", 100),
            (TssStrand::Minus, "gene-b", "SYN?A", "promoter-minus", 200)].into_iter().map(|(strand, gene, symbol, promoter, start)| {
            let tracks = matrices.iter().zip(&signals).map(|(matrix, (forward, reverse))| {
                let forward_peaks = peaks(forward);
                let reverse_peaks = peaks(reverse);
                TssProfileTrack {
                    accession: matrix.specification.source_id.clone(), motif_length_bp: matrix.matrix_counts.len(),
                    forward_scores: forward.clone(), reverse_scores: reverse.clone(),
                    forward_maximum: forward_peaks.first().cloned(), reverse_maximum: reverse_peaks.first().cloned(),
                    forward_peaks, reverse_peaks, normalization_reference: json!({"origin": "synthetic; not scored"}),
                }
            }).collect();
            let comparisons = [(TssStrand::Plus, 1.0), (TssStrand::Minus, -1.0)].into_iter().map(|(local_strand, correlation)| TssMatrixComparison {
                factor_id: "SYNTH_A".into(), left_accession: "MA0001.1".into(), right_accession: "MA0001.2".into(),
                local_strand, paired_window_count: 2, excluded_window_count: 2,
                pearson: Some(correlation), spearman: Some(correlation), undefined_reason: None,
                method: "Pearson and Spearman; same local strand; raw unclipped scores; no smoothing; common valid starts".into(),
            }).collect();
            TssProfileWindow {
                detail_context: None,
                record: TssRecord {
                    promoter_id: promoter.into(), gene_id: gene.into(), gene_symbol: symbol.into(),
                    geometry: TssGeometry { chromosome: "synthetic-contig".into(), strand, tss_1based: start + 2,
                        start_1based: start, end_1based: start + 4, upstream_bp: 2, downstream_bp: 2 },
                    transcripts: vec![format!("{promoter}-tx1"), format!("{promoter}-tx2")],
                    sequence_sha256: sha256_hex_bytes(b"ACGTA"),
                }, selected: strand == TssStrand::Plus, selection_evidence: None, tracks, comparisons,
            }
        }).collect();
        TssProfileReport {
            schema: REPORT_SCHEMA.into(),
            reference: TssReference {
                genome_id: "synthetic-genome".into(),
                assembly: "synthetic-assembly-v1".into(),
                annotation_release: Some("synthetic-annotation-v1".into()),
            },
            panel_resolution: TssPanelResolution {
                panel,
                panel_sha256,
                registry_sources: vec![binding(
                    "active_registry",
                    "synthetic-registry.json",
                    b"hand-crafted registry",
                )],
                registry_source_url: None,
                matrices,
            },
            inputs,
            windows,
            source: None,
            producer_revision: "synthetic-producer-revision".into(),
            producer_executable_sha256: None,
            lockfile_sha256: sha256_hex_bytes(
                b"synthetic producer lockfile, deliberately distinct from exporter",
            ),
            score_policy: BTreeMap::from([
                (
                    "origin".into(),
                    "hand-crafted fixture; not shared-scorer output".into(),
                ),
                (
                    "correlation_signal".into(),
                    "raw unclipped; no smoothing; common valid window starts".into(),
                ),
            ]),
            verification: "synthetic_internal_consistency_only; reference_not_assessed".into(),
            warnings: vec!["Export-only synthetic fixture".into()],
            non_claims: NON_CLAIMS.into(),
        }
    }

    fn temporary_root() -> (tempfile::TempDir, PathBuf) {
        let temp = tempfile::tempdir().unwrap();
        // macOS's system /var or /tmp aliases are not test output symlinks.
        let root = fs::canonicalize(temp.path()).unwrap();
        (temp, root)
    }

    fn synthetic_target_report() -> TssProfileReport {
        let mut report = synthetic_report();
        report.inputs[0].role = "bundle_manifest".into();
        report.source = Some(TssBundleSource {
            schema: TARGET_BUNDLE_SCHEMA.into(),
            manifest_sha256: manifest_digest(&report.inputs).unwrap().into(),
            source_revision: Some("synthetic-bundle-revision".into()),
            dataset_id: Some("synthetic-target-dataset".into()),
            producer_sha256: Some(sha256_hex_bytes(b"synthetic bundle producer bytes")),
        });
        report.reference.annotation_release = None;
        report.producer_executable_sha256 =
            Some(sha256_hex_bytes(b"synthetic profile producer bytes"));
        report.windows[0].selection_evidence = Some(TssSelectionEvidence {
            label: "Selected in synthetic report".into(),
            legend: "Synthetic descriptive window criterion; not evidence of TSS usage, direct binding or promoter activity.".into(),
            criterion: "Synthetic window overlap rule (fixture, not an assay)".into(),
            factor: Some("SYNTH_A".into()),
        });
        report
    }

    fn request(root: &Path, name: &str) -> ExportTssProfilesRequest {
        ExportTssProfilesRequest {
            context_manifest: None,
            output_dir: root.join(name).to_str().unwrap().into(),
            rendering: TssProfileRenderOptions::default(),
            formats: vec![TssExportFormat::Svg],
        }
    }

    fn files(root: &Path) -> Vec<String> {
        let mut files = fs::read_dir(root)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_str().unwrap().to_owned())
            .collect::<Vec<_>>();
        files.sort();
        files
    }

    fn score_row(
        text: &str,
        accession: &str,
        start: &str,
        strand: &str,
    ) -> BTreeMap<String, String> {
        let mut lines = text.lines().filter(|line| !line.starts_with('#'));
        let columns: Vec<_> = lines.next().unwrap().split('\t').collect();
        lines
            .map(|line| {
                columns
                    .iter()
                    .copied()
                    .zip(line.split('\t'))
                    .map(|(k, v)| (k.to_owned(), v.to_owned()))
                    .collect::<BTreeMap<_, _>>()
            })
            .find(|row| {
                row["accession"] == accession
                    && row["local_window_start_0based"] == start
                    && row["local_motif_strand"] == strand
            })
            .unwrap()
    }

    #[test]
    fn synthetic_export_preserves_subsets_raw_display_coordinates_and_inventory() {
        let (_temp, root) = temporary_root();
        let report = synthetic_report();
        let request = request(&root, "export");
        let receipt = export_tss_profiles(&report, &request).unwrap();
        let output = Path::new(&request.output_dir);
        let index: ExportIndex = read_json(&output.join(INDEX_FILE), MAX_METADATA_BYTES).unwrap();
        assert_eq!(receipt.tss_count, 2);
        assert_eq!(receipt.page_count, 2);
        assert_eq!(receipt.outputs.len(), 13);
        assert_eq!(files(output).len(), receipt.outputs.len() + 1);
        assert!(!receipt.outputs.contains_key(RECEIPT_FILE));
        assert!(!index.outputs.contains_key(INDEX_FILE));
        assert!(!index.outputs.contains_key(RECEIPT_FILE));
        assert_eq!(receipt.input_manifest_sha256, report.inputs[0].sha256);
        assert!(receipt.source_revision.is_none());
        assert!(index.source.is_none());
        assert!(index.source_revision.is_none());
        assert!(index.producer_executable_sha256.is_none());
        assert!(
            !receipt
                .inputs
                .iter()
                .any(|binding| binding.role == "producer_executable")
        );
        assert_eq!(
            receipt.lockfile_sha256,
            sha256_hex_bytes(LOCKFILE.as_bytes())
        );
        assert_eq!(index.producer_lockfile_sha256, report.lockfile_sha256);
        assert_ne!(receipt.lockfile_sha256, index.producer_lockfile_sha256);
        assert_eq!(
            receipt.exporter_revision,
            option_env!("GENTLE_SOURCE_REVISION").unwrap_or("unknown")
        );
        assert_eq!(
            receipt.executable_sha256,
            sha256_file_hex(&std::env::current_exe().unwrap()).unwrap()
        );
        assert!(
            receipt
                .inputs
                .iter()
                .any(|b| b.role == "normalized_sequence")
        );
        assert!(receipt.inputs.iter().any(|b| b.role == "resolved_matrix"));
        assert!(receipt.inputs.iter().any(|b| b.role == "active_registry"));
        for gene in &index.genes {
            assert!(portable_name(&gene.report));
            let subset: TssProfileReport =
                read_json(&output.join(&gene.report), MAX_REPORT_BYTES).unwrap();
            assert_eq!(subset.windows.len(), 1);
            assert_eq!(subset.windows[0].record.gene_id, gene.gene_id);
            assert_eq!(subset.windows[0].comparisons.len(), 2);
            assert_eq!(subset.windows[0].record.transcripts.len(), 2);
            assert_eq!(subset.panel_resolution.matrices.len(), 3);
            assert!(subset.source.is_none());
            assert!(subset.producer_executable_sha256.is_none());
            assert!(subset.windows[0].selection_evidence.is_none());
            assert!(
                fs::read_to_string(output.join(&gene.comparisons))
                    .unwrap()
                    .contains("raw_unclipped")
            );
        }
        let plus = fs::read_to_string(output.join(&index.genes[0].scores)).unwrap();
        assert!(plus.starts_with(&format!("# non_claims: {NON_CLAIMS}\n")));
        assert_eq!(
            plus.lines().filter(|line| !line.starts_with('#')).count(),
            1 + 5 * 3 * 2
        );
        let first = score_row(&plus, "MA0001.1", "0", "+");
        for (key, value) in [
            ("tss_relative_window_start_bp", "-2"),
            ("genomic_window_start_1based", "100"),
            ("motif_length_bp", "3"),
            ("motif_genomic_start_1based", "100"),
            ("motif_genomic_end_1based", "102"),
            ("motif_5prime_genomic_1based", "100"),
            ("raw_score", "-2"),
            ("display_score", "0"),
            ("selection_label", GENERIC_SELECTION_LABEL),
            ("selection_legend", GENERIC_SELECTION_LEGEND),
            ("selection_criterion", "null"),
            ("selection_factor", "null"),
        ] {
            assert_eq!(first[key], value, "{key}");
        }
        let reverse = score_row(&plus, "MA0001.1", "0", "-");
        assert_eq!(reverse["genomic_window_start_1based"], "100");
        assert_eq!(reverse["motif_5prime_genomic_1based"], "102");
        assert_eq!(reverse["genomic_motif_strand"], "-");
        let tss = score_row(&plus, "MA0001.1", "2", "+");
        assert_eq!(tss["tss_relative_window_start_bp"], "0");
        assert_eq!(tss["genomic_window_start_1based"], "102");
        assert_eq!(tss["motif_genomic_end_1based"], "104");
        let missing = score_row(&plus, "MA0001.1", "1", "+");
        assert_eq!(missing["availability"], "unavailable");
        assert_eq!(missing["raw_score"], "null");
        assert_eq!(missing["display_score"], "null");
        assert_eq!(missing["motif_genomic_start_1based"], "101");
        let tail = score_row(&plus, "MA0001.1", "4", "+");
        assert_eq!(tail["availability"], "outside_window");
        assert_eq!(tail["genomic_window_start_1based"], "104");
        for key in [
            "motif_genomic_start_1based",
            "motif_genomic_end_1based",
            "motif_5prime_genomic_1based",
            "raw_score",
            "display_score",
        ] {
            assert_eq!(tail[key], "null");
        }
        let minus = fs::read_to_string(output.join(&index.genes[1].scores)).unwrap();
        let forward = score_row(&minus, "MA0001.1", "0", "+");
        let reverse = score_row(&minus, "MA0001.1", "0", "-");
        assert_eq!(forward["genomic_window_start_1based"], "204");
        assert_eq!(forward["motif_genomic_start_1based"], "202");
        assert_eq!(forward["motif_genomic_end_1based"], "204");
        assert_eq!(forward["motif_5prime_genomic_1based"], "204");
        assert_eq!(forward["genomic_motif_strand"], "-");
        assert_eq!(reverse["genomic_window_start_1based"], "204");
        assert_eq!(reverse["motif_5prime_genomic_1based"], "202");
        assert_eq!(reverse["genomic_motif_strand"], "+");
        for field in [
            "selection_label",
            "selection_legend",
            "selection_criterion",
            "selection_factor",
        ] {
            assert_eq!(forward[field], "null");
        }
        let minus_last = score_row(&minus, "MA0001.2", "3", "+");
        assert_eq!(minus_last["genomic_window_start_1based"], "201");
        assert_eq!(minus_last["motif_genomic_start_1based"], "200");
        let readme = fs::read_to_string(output.join(README_FILE)).unwrap();
        assert!(readme.contains(NON_CLAIMS));
        assert!(readme.contains("raster-backed"));
        assert!(readme.contains("NOT verified"));
        assert!(readme.contains("Input source revision: not_assessed"));
        assert!(readme.contains("Profile producer executable SHA-256: not_assessed"));
        assert!(readme.contains(GENERIC_SELECTION_LEGEND));
        assert_eq!(index.genes[0].selections[0].label, GENERIC_SELECTION_LABEL);
        assert!(index.genes[0].selections[0].criterion.is_none());
        assert!(index.genes[1].selections.is_empty());
        assert_eq!(
            fs::read_to_string(output.join(COMPARISON_FILE))
                .unwrap()
                .lines()
                .filter(|line| !line.starts_with('#'))
                .count(),
            5
        );
        read_and_verify_tss_profile_receipt(output).unwrap();
        assert_eq!(files(&root), vec!["export"]);
    }

    #[test]
    fn source_manifest_selection_and_distinct_producer_bindings_survive_export() {
        let (_temp, root) = temporary_root();
        let mut report = synthetic_target_report();
        // Hand-crafted source-identity stub for byte binding only, not a complete
        // target-reader acceptance fixture. Recreated here without external data.
        let manifest = json!({
            "schema": TARGET_BUNDLE_SCHEMA,
            "source_revision": report.source.as_ref().unwrap().source_revision,
            "dataset_id": report.source.as_ref().unwrap().dataset_id,
            "producer_sha256": report.source.as_ref().unwrap().producer_sha256,
        });
        let manifest_path = root.join("manifest.json");
        fs::write(
            &manifest_path,
            serde_json::to_vec_pretty(&manifest).unwrap(),
        )
        .unwrap();
        let manifest_hash = sha256_file_hex(&manifest_path).unwrap();
        report.inputs[0].sha256 = manifest_hash.clone();
        report.source.as_mut().unwrap().manifest_sha256 = manifest_hash.clone();
        let req = request(&root, "export");
        let receipt = export_tss_profiles(&report, &req).unwrap();
        let output = Path::new(&req.output_dir);
        let index: ExportIndex = read_json(&output.join(INDEX_FILE), MAX_METADATA_BYTES).unwrap();
        assert_eq!(receipt.input_manifest_sha256, manifest_hash);
        assert_eq!(
            receipt.source_revision.as_deref(),
            manifest["source_revision"].as_str()
        );
        assert_ne!(
            receipt.source_revision.as_deref(),
            Some(receipt.producer_revision.as_str())
        );
        assert_eq!(index.input_manifest_sha256, receipt.input_manifest_sha256);
        assert_eq!(index.source_revision, receipt.source_revision);
        assert!(same_json(&index.source, &report.source).unwrap());
        assert!(index.reference.annotation_release.is_none());
        assert_eq!(
            index.producer_executable_sha256,
            report.producer_executable_sha256
        );
        for (role, hash) in [
            (
                "producer_executable",
                report.producer_executable_sha256.as_ref().unwrap(),
            ),
            (
                "bundle_source_producer",
                report
                    .source
                    .as_ref()
                    .unwrap()
                    .producer_sha256
                    .as_ref()
                    .unwrap(),
            ),
        ] {
            assert!(
                receipt
                    .inputs
                    .iter()
                    .any(|binding| binding.role == role && &binding.sha256 == hash)
            );
            assert_ne!(&receipt.executable_sha256, hash);
        }
        assert_eq!(
            fs::read(output.join(REPORT_FILE)).unwrap(),
            serde_json::to_vec(&report).unwrap()
        );
        let restored: TssProfileReport =
            read_json(&output.join(REPORT_FILE), MAX_REPORT_BYTES).unwrap();
        assert!(same_json(&restored, &report).unwrap());
        let gene = &index.genes[0];
        let subset: TssProfileReport =
            read_json(&output.join(&gene.report), MAX_REPORT_BYTES).unwrap();
        assert!(same_json(&subset.source, &report.source).unwrap());
        assert_eq!(
            subset.producer_executable_sha256,
            report.producer_executable_sha256
        );
        assert_eq!(subset.windows.len(), 1);
        assert!(same_json(&subset.windows[0], &report.windows[0]).unwrap());
        let evidence = report.windows[0].selection_evidence.as_ref().unwrap();
        assert_eq!(gene.selections[0].legend, evidence.legend);
        assert_eq!(
            gene.selections[0].criterion.as_deref(),
            Some(evidence.criterion.as_str())
        );
        let scores = fs::read_to_string(output.join(&gene.scores)).unwrap();
        let row = score_row(&scores, "MA0001.1", "0", "+");
        for (name, value) in [
            ("selection_label", evidence.label.as_str()),
            ("selection_legend", evidence.legend.as_str()),
            ("selection_criterion", evidence.criterion.as_str()),
            ("selection_factor", evidence.factor.as_deref().unwrap()),
        ] {
            assert_eq!(row[name], value);
        }
        for name in [
            gene.scores.as_str(),
            gene.comparisons.as_str(),
            COMPARISON_FILE,
        ] {
            let table = fs::read_to_string(output.join(name)).unwrap();
            assert!(table.starts_with(&format!("# non_claims: {NON_CLAIMS}\n")));
            assert!(table.contains("# annotation_release: null\n"));
            assert!(table.contains(&format!("# selection_legend: {}\n", evidence.legend)));
            assert!(table.contains(&format!("# input_manifest_sha256: {manifest_hash}\n")));
            assert!(table.contains(
                "selection_label\tselection_legend\tselection_criterion\tselection_factor"
            ));
        }
        let methods = fs::read_to_string(output.join(README_FILE)).unwrap();
        assert!(methods.contains(&format!("Input manifest SHA-256: {manifest_hash}")));
        assert!(methods.contains("Separately declared annotation release: not_assessed"));
        assert!(methods.contains(&evidence.label));
        assert!(methods.contains(&evidence.legend));
        assert!(methods.contains(&evidence.criterion));
        assert!(methods.contains("No CUT&RUN support is inferred."));
        verify_tss_profile_receipt(output, &receipt).unwrap();
        let replay = export_tss_profiles(&restored, &request(&root, "replay")).unwrap();
        assert_eq!(receipt.outputs, replay.outputs);
        assert!(same_json(&receipt, &replay).unwrap());
    }

    #[test]
    fn legacy_and_generic_selected_reports_do_not_invent_optional_provenance() {
        let mut report = synthetic_report();
        for role in ["manifest", "bundle_manifest"] {
            report.inputs[0].role = role.into();
            validate_tss_profile_report(&report).unwrap();
        }
        report.source = Some(TssBundleSource {
            schema: BUNDLE_SCHEMA.into(),
            manifest_sha256: manifest_digest(&report.inputs).unwrap().into(),
            source_revision: None,
            dataset_id: None,
            producer_sha256: None,
        });
        validate_tss_profile_report(&report).unwrap();
        assert!(source_revision(&report).is_none());
        let mut target = synthetic_target_report();
        target.windows[0].selection_evidence = None;
        target.producer_executable_sha256 = None;
        target.source.as_mut().unwrap().producer_sha256 = None;
        validate_tss_profile_report(&target).unwrap();
        let inputs = input_bindings(&target).unwrap();
        assert!(!inputs.iter().any(|binding| matches!(
            binding.role.as_str(),
            "producer_executable" | "bundle_source_producer"
        )));
        let methods = readme(&target).unwrap();
        assert!(methods.contains(GENERIC_SELECTION_LEGEND));
        assert!(methods.contains("Recorded criterion: not_assessed"));
    }

    #[test]
    fn malformed_source_and_selection_metadata_fail_before_writes() {
        let (_temp, root) = temporary_root();
        let mutations: Vec<Box<dyn Fn(&mut TssProfileReport)>> = vec![
            Box::new(|r| r.source.as_mut().unwrap().schema = "unrecognized".into()),
            Box::new(|r| {
                r.source.as_mut().unwrap().manifest_sha256 = sha256_hex_bytes(b"different manifest")
            }),
            Box::new(|r| r.source.as_mut().unwrap().manifest_sha256 = "A".repeat(64)),
            Box::new(|r| r.source.as_mut().unwrap().source_revision = None),
            Box::new(|r| r.source.as_mut().unwrap().source_revision = Some(" ".into())),
            Box::new(|r| {
                r.source.as_mut().unwrap().source_revision = Some("revision\nspoof".into())
            }),
            Box::new(|r| r.source.as_mut().unwrap().dataset_id = None),
            Box::new(|r| r.source.as_mut().unwrap().dataset_id = Some(" ".into())),
            Box::new(|r| r.source.as_mut().unwrap().producer_sha256 = Some("not-assessed".into())),
            Box::new(|r| r.producer_executable_sha256 = Some("A".repeat(64))),
            Box::new(|r| r.reference.annotation_release = Some("".into())),
            Box::new(|r| r.inputs[0].role = "file_with_manifest_in_its_name".into()),
            Box::new(|r| {
                r.inputs.push(TssInputBinding {
                    role: "manifest".into(),
                    name: "other-manifest.json".into(),
                    sha256: sha256_hex_bytes(b"conflicting manifest"),
                })
            }),
            Box::new(|r| r.windows[0].selected = false),
            Box::new(|r| {
                r.windows[0]
                    .selection_evidence
                    .as_mut()
                    .unwrap()
                    .label
                    .clear()
            }),
            Box::new(|r| {
                r.windows[0]
                    .selection_evidence
                    .as_mut()
                    .unwrap()
                    .legend
                    .clear()
            }),
            Box::new(|r| r.windows[0].selection_evidence.as_mut().unwrap().criterion = " ".into()),
            Box::new(|r| {
                r.windows[0].selection_evidence.as_mut().unwrap().factor = Some(" ".into())
            }),
            Box::new(|r| {
                r.windows[0].selection_evidence.as_mut().unwrap().legend =
                    "x".repeat(MAX_TEXT_BYTES + 1)
            }),
            Box::new(|r| {
                r.windows[0].selection_evidence.as_mut().unwrap().label = "invalid\0label".into()
            }),
        ];
        for (index, mutate) in mutations.into_iter().enumerate() {
            let mut report = synthetic_target_report();
            mutate(&mut report);
            assert!(
                export_tss_profiles(&report, &request(&root, "rejected")).is_err(),
                "metadata mutation {index}"
            );
            assert!(
                files(&root).is_empty(),
                "metadata mutation {index} left files"
            );
        }
        let mut legacy = synthetic_report();
        legacy.inputs[0].role = "not_a_manifest".into();
        assert!(validate_tss_profile_report(&legacy).is_err());
    }

    // Recompute hashes as a tamperer could, so semantic-linkage tests do not
    // pass merely because the ordinary file digest changed.
    fn rebind_test_output(
        output: &Path,
        receipt: &TssProfileReceipt,
        original_index: &[u8],
        name: &str,
        bytes: &[u8],
    ) -> TssProfileReceipt {
        fs::write(output.join(name), bytes).unwrap();
        let hash = sha256_hex_bytes(bytes);
        let mut altered = receipt.clone();
        altered.outputs.insert(name.into(), hash.clone());
        if name != INDEX_FILE {
            let mut index: ExportIndex = serde_json::from_slice(original_index).unwrap();
            index.outputs.insert(name.into(), hash);
            let index_bytes = serde_json::to_vec(&index).unwrap();
            fs::write(output.join(INDEX_FILE), &index_bytes).unwrap();
            altered
                .outputs
                .insert(INDEX_FILE.into(), sha256_hex_bytes(&index_bytes));
        }
        fs::write(
            output.join(RECEIPT_FILE),
            serde_json::to_vec(&altered).unwrap(),
        )
        .unwrap();
        altered
    }

    #[test]
    fn source_receipt_index_and_subset_tamper_fail_even_with_rebound_output_hashes() {
        let (_temp, root) = temporary_root();
        let report = synthetic_target_report();
        let req = request(&root, "export");
        let receipt = export_tss_profiles(&report, &req).unwrap();
        let output = Path::new(&req.output_dir);
        let receipt_bytes = fs::read(output.join(RECEIPT_FILE)).unwrap();
        let index_bytes = fs::read(output.join(INDEX_FILE)).unwrap();
        let index: ExportIndex = serde_json::from_slice(&index_bytes).unwrap();
        let subset_name = &index.genes[0].report;
        let subset_bytes = fs::read(output.join(subset_name)).unwrap();
        let mut altered = receipt.clone();
        altered.input_manifest_sha256 = sha256_hex_bytes(b"forged input manifest");
        assert!(validate_tss_profile_receipt(&altered).is_err());
        let mutations: Vec<Box<dyn Fn(&mut TssProfileReceipt)>> = vec![
            Box::new(|r| r.source_revision = Some("forged-source-revision".into())),
            Box::new(|r| r.source_revision = None),
            Box::new(|r| {
                r.input_manifest_sha256 = sha256_hex_bytes(b"forged input manifest");
                for binding in &mut r.inputs {
                    if binding.role == "bundle_manifest" {
                        binding.sha256 = r.input_manifest_sha256.clone();
                    }
                }
            }),
            Box::new(|r| {
                r.inputs
                    .iter_mut()
                    .find(|binding| binding.role == "producer_executable")
                    .unwrap()
                    .sha256 = sha256_hex_bytes(b"forged producer")
            }),
        ];
        for mutate in mutations {
            let mut altered = receipt.clone();
            mutate(&mut altered);
            validate_tss_profile_receipt(&altered).unwrap();
            fs::write(
                output.join(RECEIPT_FILE),
                serde_json::to_vec(&altered).unwrap(),
            )
            .unwrap();
            assert!(
                verify_tss_profile_receipt(output, &altered)
                    .unwrap_err()
                    .message
                    .contains("source provenance")
            );
        }
        let mutations: Vec<Box<dyn Fn(&mut ExportIndex)>> = vec![
            Box::new(|index| {
                index.input_manifest_sha256 = sha256_hex_bytes(b"forged index manifest")
            }),
            Box::new(|index| index.source_revision = Some("forged-source-revision".into())),
            Box::new(|index| {
                index.source.as_mut().unwrap().dataset_id = Some("forged-dataset".into())
            }),
            Box::new(|index| index.producer_executable_sha256 = None),
            Box::new(|index| index.reference.annotation_release = Some("inferred-release".into())),
            Box::new(|index| {
                index.genes[0].selections[0].legend = "forged selection legend".into()
            }),
        ];
        for mutate in mutations {
            let mut index: ExportIndex = serde_json::from_slice(&index_bytes).unwrap();
            mutate(&mut index);
            let altered = rebind_test_output(
                output,
                &receipt,
                &index_bytes,
                INDEX_FILE,
                &serde_json::to_vec(&index).unwrap(),
            );
            assert!(
                verify_tss_profile_receipt(output, &altered)
                    .unwrap_err()
                    .message
                    .contains("index")
            );
        }
        let mutations: Vec<Box<dyn Fn(&mut TssProfileReport)>> = vec![
            Box::new(|subset| {
                subset.source.as_mut().unwrap().source_revision =
                    Some("forged-source-revision".into())
            }),
            Box::new(|subset| subset.producer_executable_sha256 = None),
            Box::new(|subset| {
                subset.windows[0]
                    .selection_evidence
                    .as_mut()
                    .unwrap()
                    .legend = "forged legend".into()
            }),
        ];
        for mutate in mutations {
            let mut subset: TssProfileReport = serde_json::from_slice(&subset_bytes).unwrap();
            mutate(&mut subset);
            let altered = rebind_test_output(
                output,
                &receipt,
                &index_bytes,
                subset_name,
                &serde_json::to_vec(&subset).unwrap(),
            );
            assert!(
                verify_tss_profile_receipt(output, &altered)
                    .unwrap_err()
                    .message
                    .contains("unchanged subset")
            );
        }
        fs::write(output.join(subset_name), subset_bytes).unwrap();
        fs::write(output.join(INDEX_FILE), index_bytes).unwrap();
        fs::write(output.join(RECEIPT_FILE), receipt_bytes).unwrap();
        verify_tss_profile_receipt(output, &receipt).unwrap();
    }

    #[test]
    fn receipts_detect_tamper_missing_extra_files_and_empty_inventories() {
        let (_temp, root) = temporary_root();
        let request = request(&root, "export");
        let receipt = export_tss_profiles(&synthetic_report(), &request).unwrap();
        let output = Path::new(&request.output_dir);
        let original = fs::read(output.join(README_FILE)).unwrap();
        fs::write(output.join(README_FILE), b"tampered").unwrap();
        assert!(
            verify_tss_profile_receipt(output, &receipt)
                .unwrap_err()
                .message
                .contains("SHA-256 mismatch")
        );
        fs::write(output.join(README_FILE), &original).unwrap();
        fs::write(output.join("unexpected.txt"), b"extra").unwrap();
        assert!(verify_tss_profile_receipt(output, &receipt).is_err());
        fs::remove_file(output.join("unexpected.txt")).unwrap();
        fs::remove_file(output.join(README_FILE)).unwrap();
        assert!(verify_tss_profile_receipt(output, &receipt).is_err());
        fs::write(output.join(README_FILE), original).unwrap();
        let mut altered = receipt.clone();
        altered.exporter_revision = "altered-build".into();
        assert!(verify_tss_profile_receipt(output, &altered).is_err());
        altered = receipt.clone();
        altered
            .outputs
            .insert(RECEIPT_FILE.into(), sha256_hex_bytes(b"circular"));
        assert!(validate_tss_profile_receipt(&altered).is_err());
        altered = receipt.clone();
        altered.outputs.clear();
        assert!(validate_tss_profile_receipt(&altered).is_err());
        altered = receipt.clone();
        altered
            .outputs
            .insert("../outside.txt".into(), sha256_hex_bytes(b"escape"));
        assert!(validate_tss_profile_receipt(&altered).is_err());
        verify_tss_profile_receipt(output, &receipt).unwrap();
    }

    #[test]
    fn malformed_reports_fail_without_creating_files() {
        let (_temp, root) = temporary_root();
        let mutations: Vec<Box<dyn Fn(&mut TssProfileReport)>> = vec![
            Box::new(|r| r.windows.clear()),
            Box::new(|r| r.schema = "unknown".into()),
            Box::new(|r| {
                r.non_claims = "Measured binding and experimentally established TSS".into()
            }),
            Box::new(|r| r.windows[0].record.geometry.tss_1based += 1),
            Box::new(|r| r.windows[0].record.geometry.upstream_bp = usize::MAX),
            Box::new(|r| {
                r.windows[0].tracks[0].forward_scores.pop();
            }),
            Box::new(|r| r.windows[0].tracks[0].accession = "MA0001.99".into()),
            Box::new(|r| r.windows[0].tracks[0].forward_scores[0] = Some(f64::NAN)),
            Box::new(|r| r.panel_resolution.matrices[0].matrix_counts[0][0] = f64::INFINITY),
            Box::new(|r| r.panel_resolution.matrices[0].matrix_counts[0][0] = -1.0),
            Box::new(|r| r.panel_resolution.matrices[0].matrix_sha256 = "0".repeat(64)),
            Box::new(|r| {
                r.windows[0].tracks[0]
                    .forward_maximum
                    .as_mut()
                    .unwrap()
                    .score = f64::INFINITY
            }),
            Box::new(|r| r.windows[0].comparisons[0].pearson = Some(f64::NAN)),
            Box::new(|r| r.windows[0].comparisons[0].right_accession = "MA0002.1".into()),
            Box::new(|r| r.windows[0].comparisons[0].paired_window_count = usize::MAX),
            Box::new(|r| r.windows[1].record.promoter_id = r.windows[0].record.promoter_id.clone()),
            Box::new(|r| r.inputs[0].name = "../manifest.json".into()),
            Box::new(|r| {
                r.windows[0].tracks[0].normalization_reference =
                    json!("x".repeat(MAX_TEXT_BYTES + 1))
            }),
        ];
        for (index, mutate) in mutations.into_iter().enumerate() {
            let mut report = synthetic_report();
            mutate(&mut report);
            assert!(
                export_tss_profiles(&report, &request(&root, "rejected")).is_err(),
                "mutation {index}"
            );
            assert!(files(&root).is_empty(), "mutation {index} left files");
        }
    }

    #[test]
    fn comparison_inventory_and_raw_validity_counts_are_required_without_rescoring() {
        let (_temp, root) = temporary_root();
        let mutations: Vec<Box<dyn Fn(&mut TssProfileReport)>> = vec![
            Box::new(|r| r.windows[0].comparisons.clear()),
            Box::new(|r| {
                r.windows[0].comparisons.pop();
            }),
            Box::new(|r| r.windows[0].comparisons[1] = r.windows[0].comparisons[0].clone()),
            Box::new(|r| r.windows[0].comparisons[0].paired_window_count = 3),
            Box::new(|r| r.windows[0].comparisons[0].excluded_window_count = 1),
            Box::new(|r| {
                r.windows[0].tracks[0].reverse_scores[1] = Some(0.0);
                r.windows[0].tracks[1].reverse_scores[1] = Some(0.0);
            }),
        ];
        for (index, mutate) in mutations.into_iter().enumerate() {
            let mut report = synthetic_report();
            mutate(&mut report);
            assert!(
                export_tss_profiles(&report, &request(&root, "rejected")).is_err(),
                "comparison mutation {index}"
            );
            assert!(
                files(&root).is_empty(),
                "comparison mutation {index} left files"
            );
        }
        let mut report = synthetic_report();
        let window = &mut report.windows[0];
        window.tracks[0].reverse_scores[1] = Some(0.0);
        window.tracks[1].reverse_scores[1] = Some(0.0);
        window.comparisons[1].paired_window_count = 3;
        window.comparisons[1].excluded_window_count = 1;
        // The synthetic coefficients remain copied values, not independent
        // scientific acceptance evidence. This gate checks masks and identities.
        validate_tss_profile_report(&report).unwrap();
        report.windows[0].comparisons.reverse();
        for comparison in &mut report.windows[0].comparisons {
            std::mem::swap(
                &mut comparison.left_accession,
                &mut comparison.right_accession,
            );
        }
        validate_tss_profile_report(&report).unwrap();
    }

    #[test]
    fn preflight_checks_options_and_fresh_destination_without_writes() {
        let (_temp, root) = temporary_root();
        let mut req = request(&root, "new");
        preflight_tss_export(&req).unwrap();
        assert!(files(&root).is_empty());
        req.formats.clear();
        assert!(preflight_tss_export(&req).is_err());
        req.formats = vec![TssExportFormat::Svg, TssExportFormat::Svg];
        assert!(preflight_tss_export(&req).is_err());
        req.formats = vec![TssExportFormat::Svg];
        req.rendering.panels_per_page = 0;
        assert!(preflight_tss_export(&req).is_err());
        req.rendering.panels_per_page = 33;
        assert!(preflight_tss_export(&req).is_err());
        req.rendering.panels_per_page = 1;
        req.output_dir = root.join("missing/child").to_str().unwrap().into();
        assert!(preflight_tss_export(&req).is_err());
        assert!(files(&root).is_empty());
        fs::create_dir(root.join("history")).unwrap();
        fs::write(root.join("history/original.txt"), b"keep history").unwrap();
        assert!(preflight_tss_export(&request(&root, "history")).is_err());
        assert!(export_tss_profiles(&synthetic_report(), &request(&root, "history")).is_err());
        assert_eq!(
            fs::read(root.join("history/original.txt")).unwrap(),
            b"keep history"
        );
        assert_eq!(files(&root), vec!["history"]);
    }

    #[test]
    fn cancellation_cleans_staging_even_after_a_receipt_was_written() {
        let (_temp, root) = temporary_root();
        let report = synthetic_report();
        let req = request(&root, "cancelled");
        assert!(export_tss_profiles_with_cancel(&report, &req, &mut || false).is_err());
        for after_receipt in [false, true] {
            let mut saw_staging = false;
            let mut callback = || {
                for entry in fs::read_dir(&root).unwrap() {
                    let entry = entry.unwrap();
                    if entry
                        .file_name()
                        .to_string_lossy()
                        .starts_with(".gentle-tss-export-")
                        && (!after_receipt || entry.path().join(RECEIPT_FILE).exists())
                    {
                        saw_staging = true;
                        return false;
                    }
                }
                true
            };
            let error = export_tss_profiles_with_cancel(&report, &req, &mut callback).unwrap_err();
            assert!(error.message.contains("cancelled"));
            assert!(saw_staging);
            assert!(files(&root).is_empty());
        }
    }

    #[test]
    fn renderer_failure_and_shared_scale_refusal_leave_no_files() {
        let (_temp, root) = temporary_root();
        let mut report = synthetic_report();
        let mut req = request(&root, "rejected");
        req.rendering.scale_mode = Some(TssScaleMode::Shared);
        assert!(export_tss_profiles(&report, &req).is_err());
        assert!(files(&root).is_empty());
        req.rendering.scale_mode = None;
        // Each policy fits the report's text bound, but together they exceed
        // the renderer's indivisible footer height. No staging should exist.
        report.score_policy = (0..16)
            .map(|index| {
                (
                    format!("synthetic-policy-{index}"),
                    "x".repeat(MAX_TEXT_BYTES),
                )
            })
            .collect();
        validate_tss_profile_report(&report).unwrap();
        assert!(
            export_tss_profiles(&report, &req)
                .unwrap_err()
                .message
                .contains("page rendering failed")
        );
        assert!(files(&root).is_empty());
    }

    #[test]
    fn exports_are_path_independent_and_accept_an_empty_destination() {
        let (_temp, root) = temporary_root();
        fs::create_dir(root.join("first")).unwrap();
        let report = synthetic_report();
        let first = export_tss_profiles(&report, &request(&root, "first")).unwrap();
        let second = export_tss_profiles(&report, &request(&root, "second")).unwrap();
        assert_eq!(first.outputs, second.outputs);
        assert!(same_json(&first, &second).unwrap());
        assert!(export_tss_profiles(&report, &request(&root, "first")).is_err());
        assert_eq!(files(&root), vec!["first", "second"]);
        assert_ne!(gene_stem("A/B", "same"), gene_stem("A?B", "same"));
        assert_ne!(gene_stem("Gene", "ID"), gene_stem("gene", "id"));
        assert!(portable_name(&gene_stem("../CON", "a/b\\c")));
        assert_eq!(tsv("a\tb\nc\\d\r"), "a\\tb\\nc\\\\d\\r");
        assert_eq!(tsv("#gene"), "\\#gene");
        assert!(root_dependency_version("resvg").is_some());
    }

    #[test]
    fn raster_metadata_preserves_used_font_bindings_and_rejects_partial_audits() {
        // Synthetic metadata only; this test does not claim a font was rendered.
        let identities = vec![SvgUsedFontIdentity {
            families: vec!["Synthetic Family".into()],
            post_script_name: "SyntheticFamily-Regular".into(),
            face_index: 2,
            sha256: sha256_hex_bytes(b"synthetic complete font collection bytes"),
        }];
        let metadata = raster_metadata(
            &sha256_hex_bytes(b"synthetic SVG"),
            Some("page.svg"),
            "page.png",
            "png",
            100,
            80,
            7,
            &identities,
        )
        .unwrap();
        assert_eq!(
            metadata["font_identities"],
            serde_json::to_value(&identities).unwrap()
        );
        assert_eq!(metadata["used_font_face_count"], 1);
        assert_eq!(metadata["font_face_count"], 7);
        assert_eq!(metadata["font_identity_status"], FONT_IDENTITY_STATUS);
        assert_eq!(metadata["font_digest_convention"], FONT_DIGEST_CONVENTION);
        assert_eq!(
            metadata["font_reproducibility"],
            "recorded_not_enforced_or_replayed"
        );
        validate_raster_font_metadata(&metadata).unwrap();
        let mutations: Vec<Box<dyn Fn(&mut Value)>> = vec![
            Box::new(|value| value["font_identities"] = Value::Null),
            Box::new(|value| value["font_identities"] = json!([])),
            Box::new(|value| value["font_identities"][0]["sha256"] = json!("missing")),
            Box::new(|value| value["font_identities"][0]["face_index"] = json!(-1)),
            Box::new(|value| value["font_identities"][0]["families"] = Value::Null),
            Box::new(|value| value["font_identities"][0]["post_script_name"] = Value::Null),
            Box::new(|value| value["font_identities"][0]["path"] = json!("/host/private/font.otf")),
            Box::new(|value| value["used_font_face_count"] = json!(2)),
            Box::new(|value| value["font_face_count"] = json!(0)),
            Box::new(|value| value["font_identity_status"] = json!("not_audited")),
            Box::new(|value| value["font_digest_convention"] = json!("name_only")),
            Box::new(|value| value["font_reproducibility"] = json!("cross_host_replay_verified")),
            Box::new(|value| {
                let duplicate = value["font_identities"][0].clone();
                value["font_identities"]
                    .as_array_mut()
                    .unwrap()
                    .push(duplicate);
                value["used_font_face_count"] = json!(2);
            }),
        ];
        for mutate in mutations {
            let mut altered = metadata.clone();
            mutate(&mut altered);
            assert!(validate_raster_font_metadata(&altered).is_err());
        }
        assert!(
            raster_metadata(
                &sha256_hex_bytes(b"synthetic SVG"),
                None,
                "page.png",
                "png",
                100,
                80,
                7,
                &[]
            )
            .is_err()
        );
    }

    #[cfg(unix)]
    #[test]
    fn symlink_destinations_ancestors_and_artifacts_are_refused() {
        use std::os::unix::fs::symlink;
        let (_temp, root) = temporary_root();
        fs::create_dir(root.join("target")).unwrap();
        symlink(root.join("target"), root.join("link")).unwrap();
        assert!(preflight_tss_export(&request(&root, "link")).is_err());
        assert!(preflight_tss_export(&request(&root, "link/child")).is_err());
        assert!(files(&root.join("target")).is_empty());
        let receipt = export_tss_profiles(&synthetic_report(), &request(&root, "export")).unwrap();
        let output = root.join("export");
        fs::rename(output.join(README_FILE), root.join("original-readme.md")).unwrap();
        symlink(root.join("original-readme.md"), output.join(README_FILE)).unwrap();
        assert!(verify_tss_profile_receipt(&output, &receipt).is_err());
    }

    #[test]
    fn pdf_export_records_lossless_encoding_and_verifies_finished_hashes() {
        let (_temp, root) = temporary_root();
        let output = root.join("pdf-export");
        let request = ExportTssProfilesRequest {
            context_manifest: None,
            output_dir: output.to_str().unwrap().into(),
            rendering: TssProfileRenderOptions::default(),
            formats: vec![TssExportFormat::Pdf],
        };
        let report = synthetic_report();
        let receipt = export_tss_profiles(&report, &request).unwrap();
        verify_tss_profile_receipt(&output, &receipt).unwrap();
        assert_eq!(
            fs::read(output.join(REPORT_FILE)).unwrap(),
            serde_json::to_vec(&report).unwrap()
        );
        let pdfs: Vec<_> = receipt
            .render_metadata
            .iter()
            .filter(|entry| entry["format"] == "pdf")
            .collect();
        assert_eq!(pdfs.len(), receipt.page_count);
        assert!(!pdfs.is_empty());
        for entry in pdfs {
            assert_eq!(
                entry["pdf_image_encoding"],
                "FlateDecode (lossless zlib-compressed RGB)"
            );
            let bytes = fs::read(output.join(entry["output_path"].as_str().unwrap())).unwrap();
            let filter = b"/Filter /FlateDecode";
            assert!(bytes.windows(filter.len()).any(|v| v == filter));
        }
    }

    /// Opt-in visual review artifact, using the real SVG/PNG/PDF helpers. The
    /// caller supplies a fresh absolute destination so the example is retained.
    #[test]
    #[ignore = "set GENTLE_TSS_EXPORT_EXAMPLE_DIR to a fresh directory; generates retained SVG/PNG/PDF artifacts"]
    fn write_synthetic_export_example() {
        let output_dir = std::env::var("GENTLE_TSS_EXPORT_EXAMPLE_DIR")
            .expect("set GENTLE_TSS_EXPORT_EXAMPLE_DIR");
        let request = ExportTssProfilesRequest {
            context_manifest: None,
            output_dir,
            rendering: TssProfileRenderOptions::default(),
            formats: vec![
                TssExportFormat::Svg,
                TssExportFormat::Png,
                TssExportFormat::Pdf,
            ],
        };
        let receipt = export_tss_profiles(&synthetic_report(), &request).unwrap();
        let output = Path::new(&request.output_dir);
        verify_tss_profile_receipt(output, &receipt).unwrap();
        for (name, _) in &receipt.outputs {
            let metadata = receipt
                .render_metadata
                .iter()
                .find(|v| v["output_path"] == name.as_str());
            if name.ends_with(".png") || name.ends_with(".pdf") {
                let metadata = metadata.unwrap();
                validate_raster_font_metadata(metadata).unwrap();
                let identities = metadata["font_identities"].as_array().unwrap();
                assert!(!identities.is_empty());
                for identity in identities {
                    digest_field(identity["sha256"].as_str().unwrap(), "rendered font digest")
                        .unwrap();
                    assert!(identity["face_index"].as_u64().is_some());
                    assert!(identity.get("path").is_none());
                }
            } else if name.ends_with(".svg") {
                assert!(metadata.unwrap()["font_identities"].is_null());
            }
            if name.ends_with(".png") {
                assert!(
                    fs::read(output.join(name))
                        .unwrap()
                        .starts_with(b"\x89PNG\r\n\x1a\n")
                );
            } else if name.ends_with(".pdf") {
                assert!(
                    fs::read(output.join(name))
                        .unwrap()
                        .starts_with(b"%PDF-1.4")
                );
                let metadata = receipt
                    .render_metadata
                    .iter()
                    .find(|v| v["output_path"] == name.as_str())
                    .unwrap();
                assert_eq!(
                    metadata["pdf_representation"],
                    "single-page raster-backed RGB image; not vector PDF"
                );
                assert!(metadata["font_face_count"].as_u64().unwrap() > 0);
                assert!(metadata["font_identities"].is_array());
            }
        }
    }
}
