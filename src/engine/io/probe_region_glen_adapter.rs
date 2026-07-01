//! Deterministic adapter for compact Glen-style probe-region validation inputs.
//!
//! This module is intentionally narrow: it converts a small, committed CSV that
//! documents selected rows from Glen's exploratory E-MTAB-14704 analysis into
//! GENtle's canonical four-file probe-region helper-output contract. It does
//! not parse CEL files, run R/APT, or infer vendor annotation.

use super::*;

const GLEN_REGION_TABLE_FILE: &str = "region_intensity_chrom_order.csv";
const GLEN_PROBE_TABLE_FILE: &str = "probe_intensity_chrom_order.csv";
const GLEN_MATRIX_MANIFEST_FILE: &str = "normalized_feature_matrix_manifest.json";
const GLEN_PROVENANCE_FILE: &str = "provenance.json";
const GLEN_MATRIX_MANIFEST_SCHEMA: &str = "gentle.probe_region_normalized_matrix_manifest.v1";
const GLEN_BACKEND_PROVENANCE_SCHEMA: &str = "gentle.probe_region_backend_provenance.v1";

const GLEN_REQUIRED_COLUMNS: &[&str] = &[
    "probeset_id",
    "probe_id",
    "chromosome",
    "parent_start",
    "parent_stop",
    "start",
    "stop",
    "strand",
    "transcript_cluster_id",
    "gene_symbol",
    "intensity_source",
    "x",
    "y",
];

#[derive(Clone)]
struct GlenProbeRegionRow {
    probeset_id: String,
    probe_id: String,
    chromosome: String,
    parent_start: String,
    parent_stop: String,
    start: String,
    stop: String,
    strand: String,
    transcript_cluster_id: String,
    gene_symbol: String,
    intensity_source: String,
    x: String,
    y: String,
    values: Vec<String>,
}

#[derive(Default)]
struct GlenRegionAccumulator {
    chromosome: String,
    start: usize,
    stop: usize,
    strand: String,
    probeset_id: String,
    transcript_cluster_id: String,
    gene_symbol: String,
    value_sums: Vec<f64>,
    value_counts: Vec<usize>,
    probe_count: usize,
}

impl GentleEngine {
    pub(crate) fn import_glen_probe_region_fixture(
        input_path: &Path,
        output_dir: &Path,
    ) -> Result<(), EngineError> {
        let rows_and_columns = Self::glen_probe_region_rows(input_path)?;
        let (rows, value_columns) = rows_and_columns;
        if rows.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Glen probe-region adapter input '{}' contains no data rows",
                    input_path.to_string_lossy()
                ),
                cause_chain: vec![],
            });
        }

        std::fs::create_dir_all(output_dir).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create Glen probe-region output directory '{}': {e}",
                output_dir.to_string_lossy()
            ),
            cause_chain: vec![],
        })?;

        Self::write_glen_probe_table(output_dir, &rows, &value_columns)?;
        Self::write_glen_region_table(output_dir, &rows, &value_columns)?;
        Self::write_glen_manifest(output_dir, &value_columns)?;
        Self::write_glen_provenance(output_dir, input_path, rows.len(), &value_columns)?;
        Ok(())
    }

    fn glen_probe_region_rows(
        input_path: &Path,
    ) -> Result<(Vec<GlenProbeRegionRow>, Vec<String>), EngineError> {
        let mut reader = csv::ReaderBuilder::new()
            .flexible(true)
            .trim(csv::Trim::All)
            .from_path(input_path)
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not open Glen probe-region adapter input '{}': {e}",
                    input_path.to_string_lossy()
                ),
                cause_chain: vec![],
            })?;
        let headers = reader
            .headers()
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not read Glen probe-region adapter header '{}': {e}",
                    input_path.to_string_lossy()
                ),
                cause_chain: vec![],
            })?
            .clone();
        let header_map = headers
            .iter()
            .enumerate()
            .map(|(idx, header)| (header.to_string(), idx))
            .collect::<BTreeMap<_, _>>();
        let missing = GLEN_REQUIRED_COLUMNS
            .iter()
            .filter(|column| !header_map.contains_key(**column))
            .copied()
            .collect::<Vec<_>>();
        if !missing.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Glen probe-region adapter input '{}' is missing required column(s): {}",
                    input_path.to_string_lossy(),
                    missing.join(", ")
                ),
                cause_chain: vec![],
            });
        }
        let value_columns = headers
            .iter()
            .filter(|header| {
                header.ends_with(".CEL")
                    || header.starts_with("mean_log2_")
                    || header.starts_with("sd_log2_")
                    || header.starts_with("log2FC_")
            })
            .map(str::to_string)
            .collect::<Vec<_>>();
        if value_columns.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Glen probe-region adapter input '{}' has no sample, mean_log2_*, sd_log2_*, or log2FC_* value columns",
                    input_path.to_string_lossy()
                ),
                cause_chain: vec![],
            });
        }
        let value_indices = value_columns
            .iter()
            .filter_map(|column| header_map.get(column).copied())
            .collect::<Vec<_>>();

        let mut rows = Vec::new();
        for record in reader.records() {
            let record = record.map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not read Glen probe-region adapter row '{}': {e}",
                    input_path.to_string_lossy()
                ),
                cause_chain: vec![],
            })?;
            let cell = |column: &str| {
                header_map
                    .get(column)
                    .and_then(|idx| record.get(*idx))
                    .unwrap_or("")
                    .trim()
                    .to_string()
            };
            let values = value_indices
                .iter()
                .map(|idx| record.get(*idx).unwrap_or("").trim().to_string())
                .collect::<Vec<_>>();
            rows.push(GlenProbeRegionRow {
                probeset_id: cell("probeset_id"),
                probe_id: cell("probe_id"),
                chromosome: cell("chromosome"),
                parent_start: cell("parent_start"),
                parent_stop: cell("parent_stop"),
                start: cell("start"),
                stop: cell("stop"),
                strand: cell("strand"),
                transcript_cluster_id: cell("transcript_cluster_id"),
                gene_symbol: cell("gene_symbol"),
                intensity_source: cell("intensity_source"),
                x: cell("x"),
                y: cell("y"),
                values,
            });
        }
        rows.sort_by(|left, right| {
            left.chromosome
                .cmp(&right.chromosome)
                .then(Self::glen_usize(&left.start).cmp(&Self::glen_usize(&right.start)))
                .then(left.probeset_id.cmp(&right.probeset_id))
                .then(left.probe_id.cmp(&right.probe_id))
        });
        Ok((rows, value_columns))
    }

    fn write_glen_probe_table(
        output_dir: &Path,
        rows: &[GlenProbeRegionRow],
        value_columns: &[String],
    ) -> Result<(), EngineError> {
        let path = output_dir.join(GLEN_PROBE_TABLE_FILE);
        let mut writer = csv::Writer::from_path(&path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write Glen probe table '{}': {e}", path.display()),
            cause_chain: vec![],
        })?;
        let mut header = vec![
            "chromosome",
            "start",
            "stop",
            "strand",
            "probe_id",
            "x",
            "y",
            "parent_probeset_or_region_id",
            "transcript_cluster_id",
            "gene_symbol",
            "intensity_source",
        ]
        .into_iter()
        .map(str::to_string)
        .collect::<Vec<_>>();
        header.extend(value_columns.iter().cloned());
        writer.write_record(&header).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write Glen probe header '{}': {e}",
                path.display()
            ),
            cause_chain: vec![],
        })?;
        for row in rows {
            let mut record = vec![
                row.chromosome.clone(),
                row.start.clone(),
                row.stop.clone(),
                row.strand.clone(),
                row.probe_id.clone(),
                row.x.clone(),
                row.y.clone(),
                row.probeset_id.clone(),
                row.transcript_cluster_id.clone(),
                row.gene_symbol.clone(),
                row.intensity_source.clone(),
            ];
            record.extend(row.values.iter().cloned());
            writer.write_record(&record).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write Glen probe row '{}': {e}", path.display()),
                cause_chain: vec![],
            })?;
        }
        writer.flush().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not flush Glen probe table '{}': {e}", path.display()),
            cause_chain: vec![],
        })?;
        Ok(())
    }

    fn write_glen_region_table(
        output_dir: &Path,
        rows: &[GlenProbeRegionRow],
        value_columns: &[String],
    ) -> Result<(), EngineError> {
        let path = output_dir.join(GLEN_REGION_TABLE_FILE);
        let mut groups = BTreeMap::<String, GlenRegionAccumulator>::new();
        for row in rows {
            let entry =
                groups
                    .entry(row.probeset_id.clone())
                    .or_insert_with(|| GlenRegionAccumulator {
                        chromosome: row.chromosome.clone(),
                        start: Self::glen_usize(&row.parent_start),
                        stop: Self::glen_usize(&row.parent_stop),
                        strand: row.strand.clone(),
                        probeset_id: row.probeset_id.clone(),
                        transcript_cluster_id: row.transcript_cluster_id.clone(),
                        gene_symbol: row.gene_symbol.clone(),
                        value_sums: vec![0.0; value_columns.len()],
                        value_counts: vec![0; value_columns.len()],
                        probe_count: 0,
                    });
            entry.start = entry.start.min(Self::glen_usize(&row.parent_start));
            entry.stop = entry.stop.max(Self::glen_usize(&row.parent_stop));
            entry.probe_count += 1;
            for (idx, raw) in row.values.iter().enumerate() {
                if let Ok(value) = raw.parse::<f64>()
                    && value.is_finite()
                {
                    entry.value_sums[idx] += value;
                    entry.value_counts[idx] += 1;
                }
            }
        }

        let mut writer = csv::Writer::from_path(&path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write Glen region table '{}': {e}",
                path.display()
            ),
            cause_chain: vec![],
        })?;
        let mut header = vec![
            "chromosome",
            "start",
            "stop",
            "strand",
            "probeset_or_region_id",
            "transcript_cluster_id",
            "number_of_probes",
            "gene_symbol",
        ]
        .into_iter()
        .map(str::to_string)
        .collect::<Vec<_>>();
        header.extend(value_columns.iter().cloned());
        writer.write_record(&header).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write Glen region header '{}': {e}",
                path.display()
            ),
            cause_chain: vec![],
        })?;
        for group in groups.values() {
            let mut record = vec![
                group.chromosome.clone(),
                group.start.to_string(),
                group.stop.to_string(),
                group.strand.clone(),
                group.probeset_id.clone(),
                group.transcript_cluster_id.clone(),
                group.probe_count.to_string(),
                group.gene_symbol.clone(),
            ];
            for (sum, count) in group.value_sums.iter().zip(&group.value_counts) {
                if *count == 0 {
                    record.push(String::new());
                } else {
                    record.push(format!("{:.3}", sum / *count as f64));
                }
            }
            writer.write_record(&record).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write Glen region row '{}': {e}", path.display()),
                cause_chain: vec![],
            })?;
        }
        writer.flush().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not flush Glen region table '{}': {e}",
                path.display()
            ),
            cause_chain: vec![],
        })?;
        Ok(())
    }

    fn write_glen_manifest(output_dir: &Path, value_columns: &[String]) -> Result<(), EngineError> {
        let path = output_dir.join(GLEN_MATRIX_MANIFEST_FILE);
        let sample_columns = value_columns
            .iter()
            .filter(|column| column.ends_with(".CEL"))
            .cloned()
            .collect::<Vec<_>>();
        let manifest = json!({
            "schema": GLEN_MATRIX_MANIFEST_SCHEMA,
            "platform": "Clariom_D_Human",
            "platform_package": "pd.clariom.d.human",
            "coordinate_system": "hg38",
            "genome_build": "GRCh38.p14",
            "normalization": "glen-derived-rma",
            "source_accession": "E-MTAB-14704",
            "targets": ["probeset", "pm_probe"],
            "sample_columns": sample_columns,
            "artifacts": [
                GLEN_REGION_TABLE_FILE,
                GLEN_PROBE_TABLE_FILE,
                GLEN_PROVENANCE_FILE
            ]
        });
        let mut text = serde_json::to_string_pretty(&manifest).unwrap_or_else(|_| "{}".to_string());
        text.push('\n');
        std::fs::write(&path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write Glen manifest '{}': {e}", path.display()),
            cause_chain: vec![],
        })
    }

    fn write_glen_provenance(
        output_dir: &Path,
        input_path: &Path,
        row_count: usize,
        value_columns: &[String],
    ) -> Result<(), EngineError> {
        let path = output_dir.join(GLEN_PROVENANCE_FILE);
        let input_text = input_path.to_string_lossy();
        let input_sha1 = std::fs::read(input_path)
            .ok()
            .map(|bytes| format!("{:x}", Sha1::digest(&bytes)));
        let provenance = json!({
            "schema": GLEN_BACKEND_PROVENANCE_SCHEMA,
            "backend": "glen_exploratory_analysis_adapter",
            "declared_backend": "glen_exploratory_analysis_adapter",
            "selected_backend": "glen_to_gentle_probe_region_schema",
            "backend_execution_policy": "fixture_regeneration_only_no_external_execution",
            "rendered_backend_command": "none; deterministic schema adapter only",
            "adapter_function": "GentleEngine::import_glen_probe_region_fixture",
            "source_accession": "E-MTAB-14704",
            "coordinate_system": "hg38",
            "genome_build": "GRCh38.p14",
            "normalization": "glen-derived-rma",
            "probe_intensity_source": "probe_level_input",
            "source_note": "Derived from Glen's exploratory E-MTAB-14704 TP73 probe-region analysis on real TP73 GRCh38.p14 coordinates. This committed fixture uses a compact documented Glen-style adapter input; no CEL files, vendor binaries, R, or APT are run during regeneration, and the full exploratory analysis snapshot is intentionally not committed.",
            "input_schema": "gentle.glen_probe_region_adapter_input.v1",
            "input_fingerprints": [
                {
                    "path": input_text.as_ref(),
                    "role": "glen_style_adapter_input",
                    "sha1": input_sha1
                }
            ],
            "tool_version_checks": [],
            "row_count": row_count,
            "value_columns": value_columns,
            "artifacts": [
                GLEN_REGION_TABLE_FILE,
                GLEN_PROBE_TABLE_FILE,
                GLEN_MATRIX_MANIFEST_FILE
            ]
        });
        let mut text =
            serde_json::to_string_pretty(&provenance).unwrap_or_else(|_| "{}".to_string());
        text.push('\n');
        std::fs::write(&path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write Glen provenance '{}': {e}", path.display()),
            cause_chain: vec![],
        })
    }

    fn glen_usize(raw: &str) -> usize {
        raw.trim()
            .parse::<f64>()
            .ok()
            .filter(|value| value.is_finite() && *value >= 0.0)
            .map(|value| value.round() as usize)
            .unwrap_or(0)
    }
}
