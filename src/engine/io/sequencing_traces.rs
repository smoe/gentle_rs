//! Raw sequencing-trace parsing and metadata-store helpers.
//!
//! This module owns ABI/AB1/SCF intake so called-base evidence can be stored
//! in project metadata without mutating construct sequences.
//!
//! Look here for:
//! - sequencing-trace store read/write helpers
//! - ABI/AB1 and SCF parsing/detection
//! - import/list/show helpers reused by shell, CLI, and future GUI review

use super::*;

#[derive(Debug, Clone, Default)]
struct ParsedSequencingTrace {
    format: SequencingTraceFormat,
    sample_name: Option<String>,
    sample_well: Option<String>,
    run_name: Option<String>,
    machine_name: Option<String>,
    machine_model: Option<String>,
    called_bases: String,
    called_base_confidence_values: Vec<u8>,
    peak_locations: Vec<u32>,
    channel_summaries: Vec<SequencingTraceChannelSummary>,
    comments_text: Option<String>,
    warnings: Vec<String>,
}

#[derive(Debug, Clone)]
struct AbifDirectoryEntry {
    tag_name: String,
    tag_number: u32,
    element_type: u16,
    element_size: u16,
    element_count: u32,
    data_size: u32,
    data_offset: u32,
    entry_offset: usize,
}

#[derive(Debug, Clone)]
struct ScfHeader {
    samples: u32,
    samples_offset: u32,
    bases: u32,
    _bases_left_clip: u32,
    _bases_right_clip: u32,
    bases_offset: u32,
    comments_size: u32,
    comments_offset: u32,
    version: String,
    sample_size: u32,
    _code_set: u32,
    private_size: u32,
    private_offset: u32,
}

impl GentleEngine {
    pub(super) fn read_sequencing_trace_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> SequencingTraceStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<SequencingTraceStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = SEQUENCING_TRACES_SCHEMA.to_string();
        }
        store
    }

    pub(super) fn read_sequencing_trace_store(&self) -> SequencingTraceStore {
        Self::read_sequencing_trace_store_from_metadata(
            self.state.metadata.get(SEQUENCING_TRACES_METADATA_KEY),
        )
    }

    pub(super) fn write_sequencing_trace_store(
        &mut self,
        mut store: SequencingTraceStore,
    ) -> Result<(), EngineError> {
        if store.traces.is_empty() {
            self.state.metadata.remove(SEQUENCING_TRACES_METADATA_KEY);
            return Ok(());
        }
        store.schema = SEQUENCING_TRACES_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize sequencing-trace metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(SEQUENCING_TRACES_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub(super) fn upsert_sequencing_trace(
        &mut self,
        record: SequencingTraceRecord,
    ) -> Result<(), EngineError> {
        let mut store = self.read_sequencing_trace_store();
        store.traces.insert(record.trace_id.clone(), record);
        self.write_sequencing_trace_store(store)
    }

    pub(super) fn normalize_sequencing_trace_id(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "trace_id cannot be empty".to_string(),
            });
        }
        let mut out = String::with_capacity(trimmed.len());
        for ch in trimmed.chars() {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
                out.push(ch);
            } else {
                out.push('_');
            }
        }
        if out.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "trace_id must contain at least one ASCII letter, digit, '-', '_' or '.'"
                    .to_string(),
            });
        }
        Ok(out)
    }

    fn unique_sequencing_trace_id(&self, base: &str) -> String {
        let normalized =
            Self::normalize_sequencing_trace_id(base).unwrap_or_else(|_| "seq_trace".to_string());
        let store = self.read_sequencing_trace_store();
        if !store.traces.contains_key(&normalized) {
            return normalized;
        }
        let mut counter = 2usize;
        loop {
            let candidate = format!("{normalized}_{counter}");
            if !store.traces.contains_key(&candidate) {
                return candidate;
            }
            counter += 1;
        }
    }

    pub fn import_sequencing_trace(
        &mut self,
        path: &str,
        trace_id: Option<&str>,
        seq_id: Option<&str>,
    ) -> Result<SequencingTraceImportReport, EngineError> {
        let normalized_seq_id = seq_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        if let Some(seq_id) = normalized_seq_id.as_deref() {
            if !self.state.sequences.contains_key(seq_id) {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Associated sequence '{}' for sequencing trace import was not found",
                        seq_id
                    ),
                });
            }
        }
        let parsed = parse_sequencing_trace_file(path)?;
        let base_trace_id = trace_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string())
            .unwrap_or_else(|| default_trace_id_base(path, parsed.sample_name.as_deref()));
        let resolved_trace_id = if trace_id.is_some() {
            Self::normalize_sequencing_trace_id(base_trace_id.as_str())?
        } else {
            self.unique_sequencing_trace_id(base_trace_id.as_str())
        };
        let existing_record = self
            .read_sequencing_trace_store()
            .traces
            .get(resolved_trace_id.as_str())
            .cloned();
        let imported_at_unix_ms = existing_record
            .as_ref()
            .filter(|existing| existing.source_path == path)
            .map(|existing| existing.imported_at_unix_ms)
            .unwrap_or_else(Self::now_unix_ms);
        let record = SequencingTraceRecord {
            schema: SEQUENCING_TRACE_RECORD_SCHEMA.to_string(),
            trace_id: resolved_trace_id.clone(),
            format: parsed.format,
            source_path: path.to_string(),
            imported_at_unix_ms,
            seq_id: normalized_seq_id,
            sample_name: parsed.sample_name,
            sample_well: parsed.sample_well,
            run_name: parsed.run_name,
            machine_name: parsed.machine_name,
            machine_model: parsed.machine_model,
            called_bases: parsed.called_bases,
            called_base_confidence_values: parsed.called_base_confidence_values,
            peak_locations: parsed.peak_locations,
            channel_summaries: parsed.channel_summaries,
            comments_text: parsed.comments_text,
        };
        self.upsert_sequencing_trace(record.clone())?;
        Ok(import_report_from_trace_record(&record, parsed.warnings))
    }

    pub fn list_sequencing_traces(&self, seq_id_filter: Option<&str>) -> Vec<SequencingTraceSummary> {
        let filter = seq_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_ascii_lowercase());
        let mut rows = self
            .read_sequencing_trace_store()
            .traces
            .values()
            .filter(|record| {
                filter.as_ref().is_none_or(|needle| {
                    record
                        .seq_id
                        .as_deref()
                        .is_some_and(|seq_id| seq_id.to_ascii_lowercase() == *needle)
                })
            })
            .map(build_sequencing_trace_summary)
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            left.seq_id
                .as_deref()
                .unwrap_or("")
                .to_ascii_lowercase()
                .cmp(&right.seq_id.as_deref().unwrap_or("").to_ascii_lowercase())
                .then(left.imported_at_unix_ms.cmp(&right.imported_at_unix_ms))
                .then(
                    left.trace_id
                        .to_ascii_lowercase()
                        .cmp(&right.trace_id.to_ascii_lowercase()),
                )
        });
        rows
    }

    pub fn get_sequencing_trace(
        &self,
        trace_id: &str,
    ) -> Result<SequencingTraceRecord, EngineError> {
        let trace_id = Self::normalize_sequencing_trace_id(trace_id)?;
        self.read_sequencing_trace_store()
            .traces
            .get(trace_id.as_str())
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequencing trace '{}' not found", trace_id),
            })
    }

    pub(crate) fn format_sequencing_trace_summary_row(row: &SequencingTraceSummary) -> String {
        format!(
            "{} format={} seq={} bases={} channels={} source={}",
            row.trace_id,
            row.format.as_str(),
            row.seq_id.as_deref().unwrap_or("-"),
            row.called_base_count,
            row.channel_count,
            row.source_path
        )
    }

    pub(crate) fn format_sequencing_trace_detail_summary(
        record: &SequencingTraceRecord,
    ) -> String {
        format!(
            "{} format={} seq={} bases={} confidences={} peaks={} channels={} source={}",
            record.trace_id,
            record.format.as_str(),
            record.seq_id.as_deref().unwrap_or("-"),
            record.called_bases.len(),
            record.called_base_confidence_values.len(),
            record.peak_locations.len(),
            record.channel_summaries.len(),
            record.source_path
        )
    }
}

fn default_trace_id_base(path: &str, sample_name: Option<&str>) -> String {
    let sample_name = sample_name
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| value.to_string());
    sample_name.unwrap_or_else(|| {
        Path::new(path)
            .file_stem()
            .and_then(|value| value.to_str())
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("seq_trace")
            .to_string()
    })
}

fn build_sequencing_trace_summary(record: &SequencingTraceRecord) -> SequencingTraceSummary {
    SequencingTraceSummary {
        trace_id: record.trace_id.clone(),
        format: record.format,
        source_path: record.source_path.clone(),
        imported_at_unix_ms: record.imported_at_unix_ms,
        seq_id: record.seq_id.clone(),
        sample_name: record.sample_name.clone(),
        run_name: record.run_name.clone(),
        called_base_count: record.called_bases.len(),
        confidence_value_count: record.called_base_confidence_values.len(),
        peak_location_count: record.peak_locations.len(),
        channel_count: record.channel_summaries.len(),
    }
}

fn import_report_from_trace_record(
    record: &SequencingTraceRecord,
    warnings: Vec<String>,
) -> SequencingTraceImportReport {
    SequencingTraceImportReport {
        schema: SEQUENCING_TRACE_IMPORT_REPORT_SCHEMA.to_string(),
        trace_id: record.trace_id.clone(),
        format: record.format,
        source_path: record.source_path.clone(),
        imported_at_unix_ms: record.imported_at_unix_ms,
        seq_id: record.seq_id.clone(),
        sample_name: record.sample_name.clone(),
        run_name: record.run_name.clone(),
        called_base_count: record.called_bases.len(),
        confidence_value_count: record.called_base_confidence_values.len(),
        peak_location_count: record.peak_locations.len(),
        channel_count: record.channel_summaries.len(),
        warnings,
    }
}

fn parse_sequencing_trace_file(path: &str) -> Result<ParsedSequencingTrace, EngineError> {
    let bytes = std::fs::read(path).map_err(|e| EngineError {
        code: ErrorCode::Io,
        message: format!("Could not read sequencing trace '{path}': {e}"),
    })?;
    if bytes.is_empty() {
        return Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Sequencing trace '{path}' is empty"),
        });
    }
    match detect_sequencing_trace_format(&bytes)? {
        SequencingTraceFormat::AbiAb1 => parse_abi_trace(path, &bytes),
        SequencingTraceFormat::Scf => parse_scf_trace(path, &bytes),
    }
}

fn detect_sequencing_trace_format(bytes: &[u8]) -> Result<SequencingTraceFormat, EngineError> {
    if bytes.len() >= 4 && &bytes[..4] == b"ABIF" {
        return Ok(SequencingTraceFormat::AbiAb1);
    }
    if bytes.len() >= 4 && &bytes[..4] == b".scf" {
        return Ok(SequencingTraceFormat::Scf);
    }
    Err(EngineError {
        code: ErrorCode::InvalidInput,
        message:
            "Unsupported sequencing trace format (expected ABIF/AB1 or SCF magic bytes)"
                .to_string(),
    })
}

fn parse_abi_trace(path: &str, bytes: &[u8]) -> Result<ParsedSequencingTrace, EngineError> {
    let entries = parse_abif_directory(bytes)?;
    let called_entry = find_preferred_abif_entry(&entries, "PBAS", &[2, 1]).ok_or_else(|| {
        EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "ABIF trace '{path}' lacks PBAS1/PBAS2 called bases and is not a supported sequencing trace"
            ),
        }
    })?;
    let called_bases = decode_abif_text(called_entry, bytes)?;
    if called_bases.trim().is_empty() {
        return Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("ABIF trace '{path}' contains an empty PBAS called-base string"),
        });
    }
    let confidence_values = find_preferred_abif_entry(&entries, "PCON", &[2, 1])
        .map(|entry| decode_abif_u8_vec(entry, bytes))
        .transpose()?
        .unwrap_or_default();
    let peak_locations = find_preferred_abif_entry(&entries, "PLOC", &[2, 1])
        .map(|entry| decode_abif_peak_locations(entry, bytes))
        .transpose()?
        .unwrap_or_default();
    let base_order = find_abif_entry(&entries, "FWO_", 1)
        .map(|entry| decode_abif_text(entry, bytes))
        .transpose()?
        .unwrap_or_else(|| "GATC".to_string());
    let channel_summaries = abi_channel_summaries(&entries, base_order.as_str());
    let sample_name = find_abif_entry(&entries, "SMPL", 1)
        .map(|entry| decode_abif_text(entry, bytes))
        .transpose()?;
    let sample_well = find_abif_entry(&entries, "TUBE", 1)
        .map(|entry| decode_abif_text(entry, bytes))
        .transpose()?;
    let run_name = find_abif_entry(&entries, "RunN", 1)
        .map(|entry| decode_abif_text(entry, bytes))
        .transpose()?;
    let machine_name = find_abif_entry(&entries, "MCHN", 1)
        .map(|entry| decode_abif_text(entry, bytes))
        .transpose()?;
    let machine_model = find_abif_entry(&entries, "MODL", 1)
        .map(|entry| decode_abif_text(entry, bytes))
        .transpose()?;
    let comments_text = find_abif_entry(&entries, "CMNT", 1)
        .map(|entry| decode_abif_text(entry, bytes))
        .transpose()?;
    let mut warnings = vec![];
    if confidence_values.is_empty() {
        warnings.push("Trace lacks PCON1/PCON2 called-base confidence values".to_string());
    } else if confidence_values.len() != called_bases.len() {
        warnings.push(format!(
            "Called-base confidence count {} does not match called-base count {}",
            confidence_values.len(),
            called_bases.len()
        ));
    }
    if peak_locations.is_empty() {
        warnings.push("Trace lacks PLOC1/PLOC2 peak locations".to_string());
    } else if peak_locations.len() != called_bases.len() {
        warnings.push(format!(
            "Peak-location count {} does not match called-base count {}",
            peak_locations.len(),
            called_bases.len()
        ));
    }
    if channel_summaries.is_empty() {
        warnings.push("Trace lacks a complete DATA1-4 or DATA9-12 channel set".to_string());
    }
    Ok(ParsedSequencingTrace {
        format: SequencingTraceFormat::AbiAb1,
        sample_name,
        sample_well,
        run_name,
        machine_name,
        machine_model,
        called_bases,
        called_base_confidence_values: confidence_values,
        peak_locations,
        channel_summaries,
        comments_text,
        warnings,
    })
}

fn parse_scf_trace(path: &str, bytes: &[u8]) -> Result<ParsedSequencingTrace, EngineError> {
    let header = parse_scf_header(bytes)?;
    validate_scf_offsets(path, bytes, &header)?;
    let (called_bases, peak_locations, confidence_values) =
        parse_scf_base_section(bytes, &header)?;
    let channel_summaries = vec![
        SequencingTraceChannelSummary {
            channel: "A".to_string(),
            trace_set: "trace".to_string(),
            point_count: header.samples as usize,
        },
        SequencingTraceChannelSummary {
            channel: "C".to_string(),
            trace_set: "trace".to_string(),
            point_count: header.samples as usize,
        },
        SequencingTraceChannelSummary {
            channel: "G".to_string(),
            trace_set: "trace".to_string(),
            point_count: header.samples as usize,
        },
        SequencingTraceChannelSummary {
            channel: "T".to_string(),
            trace_set: "trace".to_string(),
            point_count: header.samples as usize,
        },
    ];
    let comments_text = if header.comments_size > 0 {
        let start = header.comments_offset as usize;
        let end = start + header.comments_size as usize;
        Some(
            String::from_utf8_lossy(slice_range(bytes, start, end, path, "SCF comments")?)
                .trim_end_matches('\0')
                .to_string(),
        )
    } else {
        None
    };
    let comment_map = parse_scf_comment_map(comments_text.as_deref());
    let sample_name = comment_map.get("NAME").cloned();
    let run_name = comment_map
        .get("RUN")
        .cloned()
        .or_else(|| comment_map.get("RUNN").cloned());
    let machine_model = comment_map.get("MACH").cloned();
    let sample_well = comment_map.get("LANE").cloned();
    let mut warnings = vec![];
    if confidence_values.is_empty() {
        warnings.push("SCF trace did not provide base-probability confidence values".to_string());
    } else if confidence_values.len() != called_bases.len() {
        warnings.push(format!(
            "Called-base confidence count {} does not match called-base count {}",
            confidence_values.len(),
            called_bases.len()
        ));
    }
    if peak_locations.len() != called_bases.len() {
        warnings.push(format!(
            "Peak-location count {} does not match called-base count {}",
            peak_locations.len(),
            called_bases.len()
        ));
    }
    Ok(ParsedSequencingTrace {
        format: SequencingTraceFormat::Scf,
        sample_name,
        sample_well,
        run_name,
        machine_name: None,
        machine_model,
        called_bases,
        called_base_confidence_values: confidence_values,
        peak_locations,
        channel_summaries,
        comments_text,
        warnings,
    })
}

fn parse_abif_directory(bytes: &[u8]) -> Result<Vec<AbifDirectoryEntry>, EngineError> {
    if bytes.len() < 34 {
        return Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: "ABIF trace is too small to contain a root directory".to_string(),
        });
    }
    if &bytes[..4] != b"ABIF" {
        return Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: "ABIF trace does not start with the 'ABIF' magic bytes".to_string(),
        });
    }
    let root = parse_abif_directory_entry(bytes, 6)?;
    let dir_offset = if root.data_size <= 4 {
        root.entry_offset + 20
    } else {
        root.data_offset as usize
    };
    let elem_size = root.element_size as usize;
    if elem_size < 28 {
        return Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "ABIF directory entry size {} is smaller than the 28-byte minimum",
                elem_size
            ),
        });
    }
    let mut entries = Vec::with_capacity(root.element_count as usize);
    for index in 0..root.element_count as usize {
        let offset = dir_offset
            .checked_add(index.saturating_mul(elem_size))
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: "ABIF directory offset overflow".to_string(),
            })?;
        entries.push(parse_abif_directory_entry(bytes, offset)?);
    }
    Ok(entries)
}

fn parse_abif_directory_entry(
    bytes: &[u8],
    offset: usize,
) -> Result<AbifDirectoryEntry, EngineError> {
    let raw = slice_range(
        bytes,
        offset,
        offset + 28,
        "ABIF",
        "directory entry payload",
    )?;
    Ok(AbifDirectoryEntry {
        tag_name: String::from_utf8_lossy(&raw[0..4]).to_string(),
        tag_number: u32::from_be_bytes(raw[4..8].try_into().expect("length checked")),
        element_type: u16::from_be_bytes(raw[8..10].try_into().expect("length checked")),
        element_size: u16::from_be_bytes(raw[10..12].try_into().expect("length checked")),
        element_count: u32::from_be_bytes(raw[12..16].try_into().expect("length checked")),
        data_size: u32::from_be_bytes(raw[16..20].try_into().expect("length checked")),
        data_offset: u32::from_be_bytes(raw[20..24].try_into().expect("length checked")),
        entry_offset: offset,
    })
}

fn find_abif_entry<'a>(
    entries: &'a [AbifDirectoryEntry],
    tag_name: &str,
    tag_number: u32,
) -> Option<&'a AbifDirectoryEntry> {
    entries
        .iter()
        .find(|entry| entry.tag_name == tag_name && entry.tag_number == tag_number)
}

fn find_preferred_abif_entry<'a>(
    entries: &'a [AbifDirectoryEntry],
    tag_name: &str,
    preferred_numbers: &[u32],
) -> Option<&'a AbifDirectoryEntry> {
    for number in preferred_numbers {
        if let Some(entry) = find_abif_entry(entries, tag_name, *number) {
            return Some(entry);
        }
    }
    entries
        .iter()
        .filter(|entry| entry.tag_name == tag_name)
        .max_by_key(|entry| entry.tag_number)
}

fn abi_channel_summaries(
    entries: &[AbifDirectoryEntry],
    base_order: &str,
) -> Vec<SequencingTraceChannelSummary> {
    let processed = [9u32, 10, 11, 12];
    let raw = [1u32, 2, 3, 4];
    let (trace_set, numbers) = if processed
        .iter()
        .all(|number| find_abif_entry(entries, "DATA", *number).is_some())
    {
        ("processed", processed.as_slice())
    } else if raw
        .iter()
        .all(|number| find_abif_entry(entries, "DATA", *number).is_some())
    {
        ("raw", raw.as_slice())
    } else {
        return vec![];
    };
    let order = base_order.chars().collect::<Vec<_>>();
    numbers
        .iter()
        .enumerate()
        .filter_map(|(idx, number)| {
            find_abif_entry(entries, "DATA", *number).map(|entry| {
                let point_count = if entry.element_size == 0 {
                    0
                } else {
                    (entry.data_size as usize) / (entry.element_size as usize)
                };
                let channel = order
                    .get(idx)
                    .copied()
                    .map(|ch| ch.to_string())
                    .unwrap_or_else(|| format!("channel_{}", idx + 1));
                SequencingTraceChannelSummary {
                    channel,
                    trace_set: trace_set.to_string(),
                    point_count,
                }
            })
        })
        .collect()
}

fn decode_abif_text(entry: &AbifDirectoryEntry, bytes: &[u8]) -> Result<String, EngineError> {
    let raw = abif_entry_data(entry, bytes, "ABIF text tag")?;
    let text_bytes = match entry.element_type {
        18 => {
            if raw.is_empty() {
                &[][..]
            } else {
                let len = raw[0] as usize;
                let end = 1 + len.min(raw.len().saturating_sub(1));
                &raw[1..end]
            }
        }
        19 => {
            let end = raw.iter().position(|value| *value == 0).unwrap_or(raw.len());
            &raw[..end]
        }
        _ => raw,
    };
    Ok(String::from_utf8_lossy(text_bytes).trim().to_string())
}

fn decode_abif_u8_vec(entry: &AbifDirectoryEntry, bytes: &[u8]) -> Result<Vec<u8>, EngineError> {
    Ok(abif_entry_data(entry, bytes, "ABIF byte array")?.to_vec())
}

fn decode_abif_peak_locations(
    entry: &AbifDirectoryEntry,
    bytes: &[u8],
) -> Result<Vec<u32>, EngineError> {
    let raw = abif_entry_data(entry, bytes, "ABIF peak locations")?;
    match entry.element_size {
        1 => Ok(raw.iter().map(|value| *value as u32).collect()),
        2 => {
            if raw.len() % 2 != 0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "ABIF peak-location tag '{}{}' has odd byte length {} for 2-byte values",
                        entry.tag_name, entry.tag_number, raw.len()
                    ),
                });
            }
            Ok(raw
                .chunks_exact(2)
                .map(|chunk| u16::from_be_bytes([chunk[0], chunk[1]]) as u32)
                .collect())
        }
        4 => {
            if raw.len() % 4 != 0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "ABIF peak-location tag '{}{}' has byte length {} not divisible by 4",
                        entry.tag_name, entry.tag_number, raw.len()
                    ),
                });
            }
            Ok(raw
                .chunks_exact(4)
                .map(|chunk| u32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
                .collect())
        }
        other => Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "ABIF peak-location tag '{}{}' uses unsupported element size {}",
                entry.tag_name, entry.tag_number, other
            ),
        }),
    }
}

fn abif_entry_data<'a>(
    entry: &AbifDirectoryEntry,
    bytes: &'a [u8],
    label: &str,
) -> Result<&'a [u8], EngineError> {
    let offset = if entry.data_size <= 4 {
        entry.entry_offset + 20
    } else {
        entry.data_offset as usize
    };
    slice_range(
        bytes,
        offset,
        offset + entry.data_size as usize,
        "ABIF",
        label,
    )
}

fn parse_scf_header(bytes: &[u8]) -> Result<ScfHeader, EngineError> {
    if bytes.len() < 128 {
        return Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: "SCF trace is too small to contain a 128-byte header".to_string(),
        });
    }
    if &bytes[..4] != b".scf" {
        return Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: "SCF trace does not start with the '.scf' magic bytes".to_string(),
        });
    }
    let header = ScfHeader {
        samples: read_be_u32(bytes, 4, "SCF header samples")?,
        samples_offset: read_be_u32(bytes, 8, "SCF header samples_offset")?,
        bases: read_be_u32(bytes, 12, "SCF header bases")?,
        _bases_left_clip: read_be_u32(bytes, 16, "SCF header bases_left_clip")?,
        _bases_right_clip: read_be_u32(bytes, 20, "SCF header bases_right_clip")?,
        bases_offset: read_be_u32(bytes, 24, "SCF header bases_offset")?,
        comments_size: read_be_u32(bytes, 28, "SCF header comments_size")?,
        comments_offset: read_be_u32(bytes, 32, "SCF header comments_offset")?,
        version: String::from_utf8_lossy(&bytes[36..40]).to_string(),
        sample_size: read_be_u32(bytes, 40, "SCF header sample_size")?,
        _code_set: read_be_u32(bytes, 44, "SCF header code_set")?,
        private_size: read_be_u32(bytes, 48, "SCF header private_size")?,
        private_offset: read_be_u32(bytes, 52, "SCF header private_offset")?,
    };
    if header.sample_size != 1 && header.sample_size != 2 {
        return Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "SCF trace uses unsupported sample_size {} (expected 1 or 2)",
                header.sample_size
            ),
        });
    }
    Ok(header)
}

fn validate_scf_offsets(path: &str, bytes: &[u8], header: &ScfHeader) -> Result<(), EngineError> {
    let sample_bytes = scf_sample_section_len(header)?;
    let _ = slice_range(
        bytes,
        header.samples_offset as usize,
        header.samples_offset as usize + sample_bytes,
        path,
        "SCF samples section",
    )?;
    let base_bytes = scf_base_section_len(header)?;
    let _ = slice_range(
        bytes,
        header.bases_offset as usize,
        header.bases_offset as usize + base_bytes,
        path,
        "SCF base section",
    )?;
    if header.comments_size > 0 {
        let _ = slice_range(
            bytes,
            header.comments_offset as usize,
            header.comments_offset as usize + header.comments_size as usize,
            path,
            "SCF comments section",
        )?;
    }
    if header.private_size > 0 {
        let _ = slice_range(
            bytes,
            header.private_offset as usize,
            header.private_offset as usize + header.private_size as usize,
            path,
            "SCF private section",
        )?;
    }
    Ok(())
}

fn parse_scf_base_section(
    bytes: &[u8],
    header: &ScfHeader,
) -> Result<(String, Vec<u32>, Vec<u8>), EngineError> {
    let bases = header.bases as usize;
    if bases == 0 {
        return Ok((String::new(), vec![], vec![]));
    }
    let start = header.bases_offset as usize;
    if scf_is_version_3(header) {
        let peak_offset = start;
        let prob_a_offset = peak_offset + bases * 4;
        let prob_c_offset = prob_a_offset + bases;
        let prob_g_offset = prob_c_offset + bases;
        let prob_t_offset = prob_g_offset + bases;
        let bases_offset = prob_t_offset + bases;
        let spare_offset = bases_offset + bases;
        let peak_bytes = slice_range(
            bytes,
            peak_offset,
            peak_offset + bases * 4,
            "SCF",
            "v3 peak indexes",
        )?;
        let prob_a = slice_range(bytes, prob_a_offset, prob_a_offset + bases, "SCF", "prob_A")?;
        let prob_c = slice_range(bytes, prob_c_offset, prob_c_offset + bases, "SCF", "prob_C")?;
        let prob_g = slice_range(bytes, prob_g_offset, prob_g_offset + bases, "SCF", "prob_G")?;
        let prob_t = slice_range(bytes, prob_t_offset, prob_t_offset + bases, "SCF", "prob_T")?;
        let base_chars =
            slice_range(bytes, bases_offset, bases_offset + bases, "SCF", "called bases")?;
        let _ = slice_range(
            bytes,
            spare_offset,
            spare_offset + bases * 3,
            "SCF",
            "v3 spare bytes",
        )?;
        let peak_locations = peak_bytes
            .chunks_exact(4)
            .map(|chunk| u32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
            .collect::<Vec<_>>();
        let called_bases = String::from_utf8_lossy(base_chars).to_string();
        let confidence_values = (0..bases)
            .map(|idx| scf_called_base_confidence(base_chars[idx], prob_a[idx], prob_c[idx], prob_g[idx], prob_t[idx]))
            .collect::<Vec<_>>();
        Ok((called_bases, peak_locations, confidence_values))
    } else {
        let base_bytes = slice_range(bytes, start, start + bases * 12, "SCF", "base records")?;
        let mut peak_locations = Vec::with_capacity(bases);
        let mut confidence_values = Vec::with_capacity(bases);
        let mut called_bases = String::with_capacity(bases);
        for chunk in base_bytes.chunks_exact(12) {
            let peak_index = u32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
            let prob_a = chunk[4];
            let prob_c = chunk[5];
            let prob_g = chunk[6];
            let prob_t = chunk[7];
            let base = chunk[8];
            peak_locations.push(peak_index);
            confidence_values.push(scf_called_base_confidence(
                base, prob_a, prob_c, prob_g, prob_t,
            ));
            called_bases.push(base as char);
        }
        Ok((called_bases, peak_locations, confidence_values))
    }
}

fn scf_called_base_confidence(base: u8, prob_a: u8, prob_c: u8, prob_g: u8, prob_t: u8) -> u8 {
    match (base as char).to_ascii_uppercase() {
        'A' => prob_a,
        'C' => prob_c,
        'G' => prob_g,
        'T' | 'U' => prob_t,
        _ => *[prob_a, prob_c, prob_g, prob_t].iter().max().unwrap_or(&0),
    }
}

fn scf_is_version_3(header: &ScfHeader) -> bool {
    header.version.trim_start().starts_with('3')
}

fn scf_sample_section_len(header: &ScfHeader) -> Result<usize, EngineError> {
    let channels = 4usize;
    let samples = header.samples as usize;
    let sample_size = header.sample_size as usize;
    samples
        .checked_mul(channels)
        .and_then(|value| value.checked_mul(sample_size))
        .ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: "SCF sample section length overflow".to_string(),
        })
}

fn scf_base_section_len(header: &ScfHeader) -> Result<usize, EngineError> {
    let bases = header.bases as usize;
    if scf_is_version_3(header) {
        bases.checked_mul(12).ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: "SCF v3 base section length overflow".to_string(),
        })
    } else {
        bases.checked_mul(12).ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: "SCF base section length overflow".to_string(),
        })
    }
}

fn parse_scf_comment_map(text: Option<&str>) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    let Some(text) = text else {
        return out;
    };
    for raw_line in text.lines() {
        let line = raw_line.trim_end_matches('\0').trim();
        if line.is_empty() {
            continue;
        }
        let Some((key, value)) = line.split_once('=') else {
            continue;
        };
        let key = key.trim().to_ascii_uppercase();
        let value = value.trim().to_string();
        if !key.is_empty() && !value.is_empty() {
            out.insert(key, value);
        }
    }
    out
}

fn read_be_u32(bytes: &[u8], offset: usize, label: &str) -> Result<u32, EngineError> {
    let raw = slice_range(bytes, offset, offset + 4, "SCF", label)?;
    Ok(u32::from_be_bytes(
        raw.try_into().expect("length checked to 4"),
    ))
}

fn slice_range<'a>(
    bytes: &'a [u8],
    start: usize,
    end: usize,
    context: &str,
    label: &str,
) -> Result<&'a [u8], EngineError> {
    if start > end || end > bytes.len() {
        return Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "{context} {label} range {}..{} exceeds file length {}",
                start,
                end,
                bytes.len()
            ),
        });
    }
    Ok(&bytes[start..end])
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn synthetic_scf_bytes() -> Vec<u8> {
        let samples = 4u32;
        let bases = 3u32;
        let sample_size = 1u32;
        let samples_offset = 128u32;
        let sample_section_len = samples * 4 * sample_size;
        let bases_offset = samples_offset + sample_section_len;
        let base_section_len = bases * 12;
        let comments = b"NAME=synthetic_scf\nRUN=demo_run\nMACH=staden_test\n\0";
        let comments_offset = bases_offset + base_section_len;
        let total_len = comments_offset as usize + comments.len();
        let mut bytes = vec![0u8; total_len];
        bytes[0..4].copy_from_slice(b".scf");
        bytes[4..8].copy_from_slice(&samples.to_be_bytes());
        bytes[8..12].copy_from_slice(&samples_offset.to_be_bytes());
        bytes[12..16].copy_from_slice(&bases.to_be_bytes());
        bytes[16..20].copy_from_slice(&0u32.to_be_bytes());
        bytes[20..24].copy_from_slice(&0u32.to_be_bytes());
        bytes[24..28].copy_from_slice(&bases_offset.to_be_bytes());
        bytes[28..32].copy_from_slice(&(comments.len() as u32).to_be_bytes());
        bytes[32..36].copy_from_slice(&comments_offset.to_be_bytes());
        bytes[36..40].copy_from_slice(b"3.00");
        bytes[40..44].copy_from_slice(&sample_size.to_be_bytes());
        bytes[44..48].copy_from_slice(&0u32.to_be_bytes());
        bytes[48..52].copy_from_slice(&0u32.to_be_bytes());
        bytes[52..56].copy_from_slice(&0u32.to_be_bytes());

        let base_start = bases_offset as usize;
        let peak_indexes = [1u32, 2u32, 3u32];
        for (idx, peak_index) in peak_indexes.iter().enumerate() {
            let offset = base_start + idx * 4;
            bytes[offset..offset + 4].copy_from_slice(&peak_index.to_be_bytes());
        }
        let prob_a_offset = base_start + bases as usize * 4;
        let prob_c_offset = prob_a_offset + bases as usize;
        let prob_g_offset = prob_c_offset + bases as usize;
        let prob_t_offset = prob_g_offset + bases as usize;
        let called_bases_offset = prob_t_offset + bases as usize;
        bytes[prob_a_offset..prob_a_offset + 3].copy_from_slice(&[70, 1, 1]);
        bytes[prob_c_offset..prob_c_offset + 3].copy_from_slice(&[1, 71, 1]);
        bytes[prob_g_offset..prob_g_offset + 3].copy_from_slice(&[1, 1, 72]);
        bytes[prob_t_offset..prob_t_offset + 3].copy_from_slice(&[0, 0, 0]);
        bytes[called_bases_offset..called_bases_offset + 3].copy_from_slice(b"ACG");
        bytes[comments_offset as usize..comments_offset as usize + comments.len()]
            .copy_from_slice(comments);
        bytes
    }

    #[test]
    fn parse_scf_trace_extracts_called_bases_and_comment_metadata() {
        let bytes = synthetic_scf_bytes();
        let parsed = parse_scf_trace("synthetic.scf", &bytes).expect("parse SCF");
        assert_eq!(parsed.format, SequencingTraceFormat::Scf);
        assert_eq!(parsed.called_bases, "ACG");
        assert_eq!(parsed.called_base_confidence_values, vec![70, 71, 72]);
        assert_eq!(parsed.peak_locations, vec![1, 2, 3]);
        assert_eq!(parsed.sample_name.as_deref(), Some("synthetic_scf"));
        assert_eq!(parsed.run_name.as_deref(), Some("demo_run"));
        assert_eq!(parsed.machine_model.as_deref(), Some("staden_test"));
        assert_eq!(parsed.channel_summaries.len(), 4);
        assert_eq!(parsed.channel_summaries[0].point_count, 4);
    }

    #[test]
    fn import_sequencing_trace_reuses_explicit_trace_id_deterministically() {
        let td = tempdir().expect("tempdir");
        let trace_path = td.path().join("synthetic.scf");
        std::fs::write(&trace_path, synthetic_scf_bytes()).expect("write SCF");

        let mut engine = GentleEngine::default();
        let first = engine
            .import_sequencing_trace(trace_path.to_str().expect("utf-8 path"), Some("trace_a"), None)
            .expect("first import");
        let second = engine
            .import_sequencing_trace(trace_path.to_str().expect("utf-8 path"), Some("trace_a"), None)
            .expect("second import");

        assert_eq!(first.trace_id, "trace_a");
        assert_eq!(second.trace_id, "trace_a");
        assert_eq!(first.imported_at_unix_ms, second.imported_at_unix_ms);
        let stored = engine
            .get_sequencing_trace("trace_a")
            .expect("stored trace");
        assert_eq!(stored.called_bases, "ACG");
    }
}
