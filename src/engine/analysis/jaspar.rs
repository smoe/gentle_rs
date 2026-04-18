//! Shared JASPAR entry presentation and expert-inspection helpers.
//!
//! This module owns the deterministic "show me what this motif likes, what
//! random DNA does, and what JASPAR says about it" reports used by shell/CLI,
//! agents, and the dedicated GUI expert window.
//!
//! Look here for:
//! - per-entry max/min scoring sequence derivation
//! - deterministic pseudorandom background scoring summaries
//! - expert-column/logo payload derivation from local PFMs
//! - optional JASPAR REST metadata enrichment and JSON export helpers

use super::*;
use serde_json::Value;
use std::collections::BTreeMap;

impl GentleEngine {
    pub(crate) fn jaspar_remote_metadata_summary(
        metadata: &JasparRemoteMetadata,
    ) -> JasparCatalogRemoteSummary {
        let mut species_preview = metadata
            .species_assignments
            .iter()
            .map(|row| row.scientific_name.trim().to_string())
            .filter(|row| !row.is_empty())
            .collect::<Vec<_>>();
        species_preview.sort_by_key(|value| value.to_ascii_uppercase());
        species_preview.dedup_by(|left, right| left.eq_ignore_ascii_case(right));
        species_preview.truncate(3);
        JasparCatalogRemoteSummary {
            collection: metadata.collection.clone(),
            tax_group: metadata.tax_group.clone(),
            tf_class: metadata.tf_class.clone(),
            tf_family: metadata.tf_family.clone(),
            data_type: metadata.data_type.clone(),
            species_count: metadata.species_assignments.len(),
            species_preview,
        }
    }

    fn jaspar_extreme_sequence(score_matrix: &[[f64; 4]], pick_max: bool) -> String {
        const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
        let mut out = String::with_capacity(score_matrix.len());
        for column in score_matrix {
            let mut best_idx = 0usize;
            let mut best_score = column[0];
            for (idx, score) in column.iter().copied().enumerate().skip(1) {
                let better = if pick_max {
                    score > best_score
                } else {
                    score < best_score
                };
                if better {
                    best_idx = idx;
                    best_score = score;
                }
            }
            out.push(BASES[best_idx]);
        }
        out
    }

    fn deterministic_random_dna_bytes(length_bp: usize, seed: u64) -> Vec<u8> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let mut state = if seed == 0 {
            DEFAULT_JASPAR_PRESENTATION_RANDOM_SEED
        } else {
            seed
        };
        let mut out = Vec::with_capacity(length_bp);
        for _ in 0..length_bp {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            out.push(BASES[((state >> 62) & 0b11) as usize]);
        }
        out
    }

    fn scan_single_jaspar_scores(sequence: &[u8], score_matrix: &[[f64; 4]]) -> Vec<f64> {
        if score_matrix.is_empty() || sequence.len() < score_matrix.len() {
            return vec![];
        }
        let mut scores = Vec::with_capacity(
            sequence
                .len()
                .saturating_sub(score_matrix.len())
                .saturating_add(1)
                .saturating_mul(2),
        );
        let len = score_matrix.len();
        for start in 0..=(sequence.len() - len) {
            let window = &sequence[start..start + len];
            if let Some(score) = Self::score_matrix_window(window, score_matrix) {
                scores.push(score);
            }
            let rc_window = Self::reverse_complement_bytes(window);
            if let Some(score) = Self::score_matrix_window(&rc_window, score_matrix) {
                scores.push(score);
            }
        }
        scores
    }

    fn summarize_jaspar_score_distribution(scores: &[f64]) -> JasparScoreDistributionSummary {
        fn percentile(sorted_scores: &[f64], fraction: f64) -> f64 {
            if sorted_scores.is_empty() {
                return 0.0;
            }
            let capped = fraction.clamp(0.0, 1.0);
            let idx = ((sorted_scores.len().saturating_sub(1) as f64) * capped).round() as usize;
            sorted_scores[idx.min(sorted_scores.len().saturating_sub(1))]
        }

        if scores.is_empty() {
            return JasparScoreDistributionSummary::default();
        }

        let sample_count = scores.len();
        let mut sorted_scores = scores.to_vec();
        sorted_scores.sort_by(|left, right| left.total_cmp(right));

        let sum = scores.iter().sum::<f64>();
        let mean_score = sum / sample_count as f64;
        let variance = scores
            .iter()
            .map(|score| {
                let delta = *score - mean_score;
                delta * delta
            })
            .sum::<f64>()
            / sample_count as f64;

        JasparScoreDistributionSummary {
            sample_count,
            min_score: *sorted_scores.first().unwrap_or(&0.0),
            max_score: *sorted_scores.last().unwrap_or(&0.0),
            mean_score,
            stddev_score: variance.sqrt(),
            p01_score: percentile(&sorted_scores, 0.01),
            p05_score: percentile(&sorted_scores, 0.05),
            p25_score: percentile(&sorted_scores, 0.25),
            p50_score: percentile(&sorted_scores, 0.50),
            p75_score: percentile(&sorted_scores, 0.75),
            p95_score: percentile(&sorted_scores, 0.95),
            p99_score: percentile(&sorted_scores, 0.99),
            positive_fraction: scores.iter().filter(|score| **score > 0.0).count() as f64
                / sample_count as f64,
        }
    }

    fn summarize_jaspar_score_histogram(
        scores: &[f64],
        bin_count: usize,
    ) -> Vec<JasparScoreDistributionBin> {
        if scores.is_empty() || bin_count == 0 {
            return vec![];
        }
        let min_score = scores.iter().copied().fold(f64::INFINITY, f64::min);
        let max_score = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if !min_score.is_finite() || !max_score.is_finite() {
            return vec![];
        }
        if (max_score - min_score).abs() <= f64::EPSILON {
            return vec![JasparScoreDistributionBin {
                start_score: min_score,
                end_score: max_score,
                count: scores.len(),
            }];
        }

        let width = (max_score - min_score) / bin_count as f64;
        let mut counts = vec![0usize; bin_count];
        for score in scores {
            let raw_idx = ((score - min_score) / width).floor() as isize;
            let idx = raw_idx.clamp(0, bin_count.saturating_sub(1) as isize) as usize;
            counts[idx] = counts[idx].saturating_add(1);
        }
        counts
            .into_iter()
            .enumerate()
            .map(|(idx, count)| {
                let start_score = min_score + width * idx as f64;
                let end_score = if idx + 1 == bin_count {
                    max_score
                } else {
                    min_score + width * (idx + 1) as f64
                };
                JasparScoreDistributionBin {
                    start_score,
                    end_score,
                    count,
                }
            })
            .collect()
    }

    fn jaspar_score_distribution_panel(
        score_kind: TfbsScoreTrackValueKind,
        label: &str,
        score_matrix: &[[f64; 4]],
        random_background: &[u8],
    ) -> JasparScoreDistributionPanel {
        let maximizing_sequence = Self::jaspar_extreme_sequence(score_matrix, true);
        let minimizing_sequence = Self::jaspar_extreme_sequence(score_matrix, false);
        let maximizing_score =
            Self::score_matrix_window(maximizing_sequence.as_bytes(), score_matrix).unwrap_or(0.0);
        let minimizing_score =
            Self::score_matrix_window(minimizing_sequence.as_bytes(), score_matrix).unwrap_or(0.0);
        let scores = Self::scan_single_jaspar_scores(random_background, score_matrix);
        let mut sorted_scores = scores.clone();
        sorted_scores.sort_by(|left, right| left.total_cmp(right));
        JasparScoreDistributionPanel {
            score_kind,
            label: label.to_string(),
            maximizing_sequence,
            maximizing_score,
            maximizing_quantile: Self::empirical_quantile(&sorted_scores, maximizing_score),
            minimizing_sequence,
            minimizing_score,
            minimizing_quantile: Self::empirical_quantile(&sorted_scores, minimizing_score),
            distribution: Self::summarize_jaspar_score_distribution(&scores),
            histogram_bins: Self::summarize_jaspar_score_histogram(&scores, 40),
        }
    }

    fn jaspar_expert_columns(matrix_counts: &[[f64; 4]]) -> Vec<JasparExpertColumn> {
        const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
        matrix_counts
            .iter()
            .enumerate()
            .map(|(idx, counts)| {
                let total_count = counts.iter().sum::<f64>();
                let fractions = if total_count > 0.0 {
                    [
                        counts[0] / total_count,
                        counts[1] / total_count,
                        counts[2] / total_count,
                        counts[3] / total_count,
                    ]
                } else {
                    [0.0; 4]
                };
                let information_content_bits = if total_count > 0.0 {
                    let entropy = fractions
                        .iter()
                        .copied()
                        .filter(|value| *value > 0.0)
                        .map(|value| value * value.log2())
                        .sum::<f64>();
                    (2.0 + entropy).max(0.0)
                } else {
                    0.0
                };
                let dominant_idx = counts
                    .iter()
                    .copied()
                    .enumerate()
                    .max_by(|left, right| left.1.total_cmp(&right.1))
                    .map(|(idx, _)| idx)
                    .unwrap_or(0);
                JasparExpertColumn {
                    position_1based: idx + 1,
                    total_count,
                    dominant_base: BASES[dominant_idx].to_string(),
                    a_count: counts[0],
                    c_count: counts[1],
                    g_count: counts[2],
                    t_count: counts[3],
                    a_fraction: fractions[0],
                    c_fraction: fractions[1],
                    g_fraction: fractions[2],
                    t_fraction: fractions[3],
                    information_content_bits,
                    a_logo_bits: fractions[0] * information_content_bits,
                    c_logo_bits: fractions[1] * information_content_bits,
                    g_logo_bits: fractions[2] * information_content_bits,
                    t_logo_bits: fractions[3] * information_content_bits,
                }
            })
            .collect()
    }

    fn read_jaspar_text_input(path_or_url: &str) -> Result<String, EngineError> {
        if let Some(path) = path_or_url.strip_prefix("file://") {
            return std::fs::read_to_string(path).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not read JASPAR metadata file '{path}': {e}"),
            });
        }
        if path_or_url.starts_with("http://") || path_or_url.starts_with("https://") {
            let response = std::panic::catch_unwind(|| reqwest::blocking::get(path_or_url))
                .map_err(|_| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Could not fetch JASPAR metadata URL '{path_or_url}': networking backend panicked"
                    ),
                })?
                .map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Could not fetch JASPAR metadata URL '{path_or_url}': {e}"),
                })?;
            if !response.status().is_success() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Could not fetch JASPAR metadata URL '{path_or_url}': HTTP {}",
                        response.status()
                    ),
                });
            }
            return response.text().map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not read JASPAR metadata response '{path_or_url}': {e}"),
            });
        }
        std::fs::read_to_string(path_or_url).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not read JASPAR metadata file '{path_or_url}': {e}"),
        })
    }

    fn jaspar_remote_string(value: &Value) -> Option<String> {
        match value {
            Value::String(text) => {
                let trimmed = text.trim();
                (!trimmed.is_empty()).then_some(trimmed.to_string())
            }
            Value::Number(number) => Some(number.to_string()),
            Value::Object(object) => ["name", "label", "value", "scientific_name", "species"]
                .iter()
                .find_map(|key| object.get(*key).and_then(Self::jaspar_remote_string)),
            _ => None,
        }
    }

    fn jaspar_remote_lookup(value: &Value, keys: &[&str]) -> Option<String> {
        keys.iter()
            .find_map(|key| value.get(*key).and_then(Self::jaspar_remote_string))
    }

    pub(crate) fn parse_jaspar_species_assignments_value(
        value: Option<&Value>,
    ) -> Vec<JasparSpeciesAssignment> {
        let mut out = vec![];
        let Some(value) = value else {
            return out;
        };
        match value {
            Value::Array(rows) => {
                for row in rows {
                    match row {
                        Value::Object(object) => {
                            let scientific_name = ["name", "species", "scientific_name"]
                                .iter()
                                .find_map(|key| {
                                    object.get(*key).and_then(Self::jaspar_remote_string)
                                });
                            let common_name = ["common_name", "display_name", "common"]
                                .iter()
                                .find_map(|key| {
                                    object.get(*key).and_then(Self::jaspar_remote_string)
                                });
                            let tax_id = ["tax_id", "taxon_id", "id"].iter().find_map(|key| {
                                object.get(*key).and_then(Self::jaspar_remote_string)
                            });
                            if let Some(scientific_name) = scientific_name {
                                out.push(JasparSpeciesAssignment {
                                    tax_id,
                                    scientific_name,
                                    common_name,
                                });
                            }
                        }
                        Value::String(text) => {
                            let trimmed = text.trim();
                            if !trimmed.is_empty() {
                                out.push(JasparSpeciesAssignment {
                                    tax_id: None,
                                    scientific_name: trimmed.to_string(),
                                    common_name: None,
                                });
                            }
                        }
                        _ => {}
                    }
                }
            }
            Value::Object(_) | Value::String(_) => {
                let scientific_name = Self::jaspar_remote_string(value);
                if let Some(scientific_name) = scientific_name {
                    out.push(JasparSpeciesAssignment {
                        tax_id: Self::jaspar_remote_lookup(value, &["tax_id", "taxon_id", "id"]),
                        scientific_name,
                        common_name: Self::jaspar_remote_lookup(
                            value,
                            &["common_name", "display_name", "common"],
                        ),
                    });
                }
            }
            _ => {}
        }
        out.sort_by(|left, right| {
            left.scientific_name
                .to_ascii_uppercase()
                .cmp(&right.scientific_name.to_ascii_uppercase())
                .then_with(|| left.tax_id.cmp(&right.tax_id))
        });
        out.dedup_by(|left, right| {
            left.tax_id == right.tax_id
                && left
                    .scientific_name
                    .eq_ignore_ascii_case(&right.scientific_name)
        });
        out
    }

    pub(crate) fn parse_jaspar_remote_metadata_value(
        value: &Value,
        source_url: &str,
    ) -> JasparRemoteMetadata {
        JasparRemoteMetadata {
            source_url: source_url.to_string(),
            collection: Self::jaspar_remote_lookup(value, &["collection"]),
            tax_group: Self::jaspar_remote_lookup(value, &["tax_group", "taxon"]),
            tf_class: Self::jaspar_remote_lookup(value, &["tf_class", "class"]),
            tf_family: Self::jaspar_remote_lookup(value, &["tf_family", "family"]),
            data_type: Self::jaspar_remote_lookup(value, &["data_type", "experiment_type"]),
            species_assignments: Self::parse_jaspar_species_assignments_value(value.get("species")),
        }
    }

    fn jaspar_remote_metadata_url(motif_id: &str) -> String {
        format!("https://jaspar.elixir.no/api/v1/matrix/{motif_id}/?format=json")
    }

    fn fetch_jaspar_remote_metadata(
        &self,
        motif_id: &str,
    ) -> Result<JasparRemoteMetadata, EngineError> {
        let url = Self::jaspar_remote_metadata_url(motif_id);
        let text = Self::read_jaspar_text_input(&url)?;
        let value = serde_json::from_str::<Value>(&text).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not parse JASPAR metadata response for '{motif_id}': {e}"),
        })?;
        Ok(Self::parse_jaspar_remote_metadata_value(&value, &url))
    }

    fn default_jaspar_remote_metadata_snapshot_path(path: Option<&str>) -> String {
        path.map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or(crate::resource_sync::DEFAULT_JASPAR_REMOTE_METADATA_PATH)
            .to_string()
    }

    fn jaspar_remote_metadata_snapshot_row(
        motif_id: &str,
        motif_name: Option<&str>,
        consensus_iupac: &str,
        motif_length_bp: usize,
        remote_metadata: JasparRemoteMetadata,
    ) -> JasparRemoteMetadataSnapshotRow {
        JasparRemoteMetadataSnapshotRow {
            motif_id: motif_id.to_string(),
            motif_name: motif_name
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string),
            consensus_iupac: consensus_iupac.to_string(),
            motif_length_bp,
            remote_summary: Some(Self::jaspar_remote_metadata_summary(&remote_metadata)),
            remote_metadata,
        }
    }

    fn load_jaspar_remote_metadata_snapshot_rows(
        path: Option<&str>,
    ) -> Result<BTreeMap<String, JasparRemoteMetadataSnapshotRow>, String> {
        let path = Self::default_jaspar_remote_metadata_snapshot_path(path);
        let snapshot_path = std::path::Path::new(&path);
        if !snapshot_path.exists() {
            return Ok(BTreeMap::new());
        }
        let text = std::fs::read_to_string(snapshot_path)
            .map_err(|e| format!("Could not read JASPAR remote-metadata snapshot '{path}': {e}"))?;
        let snapshot =
            serde_json::from_str::<JasparRemoteMetadataSnapshot>(&text).map_err(|e| {
                format!("Could not parse JASPAR remote-metadata snapshot '{path}': {e}")
            })?;
        if !snapshot
            .schema
            .starts_with("gentle.jaspar_remote_metadata_snapshot.v")
        {
            return Err(format!(
                "Unexpected JASPAR remote-metadata snapshot schema '{}' in '{}'",
                snapshot.schema, path
            ));
        }
        Ok(snapshot
            .rows
            .into_iter()
            .map(|mut row| {
                if row.remote_summary.is_none() {
                    row.remote_summary =
                        Some(Self::jaspar_remote_metadata_summary(&row.remote_metadata));
                }
                (row.motif_id.clone(), row)
            })
            .collect())
    }

    fn select_jaspar_catalog_rows(
        &self,
        motifs: &[String],
        filter: Option<&str>,
        limit: Option<usize>,
    ) -> Result<(usize, Vec<String>, Vec<tf_motifs::TfMotifSummary>), EngineError> {
        let registry_rows = tf_motifs::list_motif_summaries();
        let registry_entry_count = registry_rows.len();
        let mut seen_ids = BTreeSet::new();
        let mut requested_motifs = vec![];
        let mut rows = vec![];
        if !motifs.is_empty() {
            for token in motifs
                .iter()
                .map(|value| value.trim())
                .filter(|value| !value.is_empty())
            {
                requested_motifs.push(token.to_string());
                let motif_definition =
                    tf_motifs::resolve_motif_definition(token).ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "JASPAR entry '{}' was not found in the local motif registry",
                            token
                        ),
                    })?;
                if !seen_ids.insert(motif_definition.id.clone()) {
                    continue;
                }
                rows.push(tf_motifs::TfMotifSummary {
                    id: motif_definition.id,
                    name: motif_definition.name,
                    consensus_iupac: motif_definition.consensus_iupac,
                    motif_length_bp: motif_definition.matrix_counts.len(),
                });
                if let Some(limit) = limit {
                    if rows.len() >= limit {
                        break;
                    }
                }
            }
            return Ok((registry_entry_count, requested_motifs, rows));
        }

        let requested_filter = filter.map(str::trim).filter(|value| !value.is_empty());
        let normalized_filter = requested_filter.map(|value| value.to_ascii_uppercase());
        let filtered_rows = registry_rows
            .into_iter()
            .filter(|row| {
                let Some(filter) = normalized_filter.as_deref() else {
                    return true;
                };
                row.id.to_ascii_uppercase().contains(filter)
                    || row
                        .name
                        .as_deref()
                        .unwrap_or("")
                        .to_ascii_uppercase()
                        .contains(filter)
                    || row.consensus_iupac.to_ascii_uppercase().contains(filter)
            })
            .take(limit.unwrap_or(usize::MAX))
            .collect::<Vec<_>>();
        Ok((registry_entry_count, requested_motifs, filtered_rows))
    }

    pub(crate) fn sync_jaspar_remote_metadata_snapshot_with_fetcher<F>(
        &self,
        motifs: &[String],
        filter: Option<&str>,
        limit: Option<usize>,
        path: Option<&str>,
        mut fetcher: F,
    ) -> Result<JasparRemoteMetadataSnapshot, EngineError>
    where
        F: FnMut(&str) -> Result<JasparRemoteMetadata, EngineError>,
    {
        let output_path = Self::default_jaspar_remote_metadata_snapshot_path(path);
        let (registry_entry_count, requested_motifs, selected_rows) =
            self.select_jaspar_catalog_rows(motifs, filter, limit)?;
        let mut warnings = vec![];
        let mut existing_rows =
            match Self::load_jaspar_remote_metadata_snapshot_rows(Some(&output_path)) {
                Ok(rows) => rows,
                Err(err) => {
                    warnings.push(err);
                    BTreeMap::new()
                }
            };
        let mut fetched_entry_count = 0usize;
        for row in selected_rows {
            match fetcher(&row.id) {
                Ok(metadata) => {
                    fetched_entry_count += 1;
                    existing_rows.insert(
                        row.id.clone(),
                        Self::jaspar_remote_metadata_snapshot_row(
                            &row.id,
                            row.name.as_deref(),
                            &row.consensus_iupac,
                            row.motif_length_bp,
                            metadata,
                        ),
                    );
                }
                Err(err) => {
                    warnings.push(format!("{}: {}", row.id, err.message));
                }
            }
        }

        let rows = existing_rows.into_values().collect::<Vec<_>>();
        let snapshot = JasparRemoteMetadataSnapshot {
            schema: JASPAR_REMOTE_METADATA_SNAPSHOT_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            filter: filter
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string),
            limit,
            requested_motifs,
            registry_entry_count,
            fetched_entry_count,
            persisted_entry_count: rows.len(),
            source: "jaspar_rest_api_v1".to_string(),
            rows,
            warnings,
        };
        Ok(snapshot)
    }

    pub(crate) fn sync_jaspar_remote_metadata_snapshot(
        &self,
        motifs: &[String],
        filter: Option<&str>,
        limit: Option<usize>,
        path: Option<&str>,
    ) -> Result<JasparRemoteMetadataSnapshot, EngineError> {
        self.sync_jaspar_remote_metadata_snapshot_with_fetcher(
            motifs,
            filter,
            limit,
            path,
            |motif_id| self.fetch_jaspar_remote_metadata(motif_id),
        )
    }

    pub(crate) fn inspect_jaspar_entry(
        &self,
        motif: &str,
        random_sequence_length_bp: usize,
        random_seed: u64,
        include_remote_metadata: bool,
        refresh_remote_metadata: bool,
    ) -> Result<JasparEntryExpertView, EngineError> {
        self.inspect_jaspar_entry_with_snapshot_path(
            motif,
            random_sequence_length_bp,
            random_seed,
            include_remote_metadata,
            refresh_remote_metadata,
            None,
        )
    }

    pub(crate) fn inspect_jaspar_entry_with_snapshot_path(
        &self,
        motif: &str,
        random_sequence_length_bp: usize,
        random_seed: u64,
        include_remote_metadata: bool,
        refresh_remote_metadata: bool,
        snapshot_path: Option<&str>,
    ) -> Result<JasparEntryExpertView, EngineError> {
        let motif = motif.trim();
        if motif.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "InspectJasparEntry requires non-empty motif".to_string(),
            });
        }
        if random_sequence_length_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "InspectJasparEntry requires random_sequence_length_bp >= 1".to_string(),
            });
        }

        let motif_definition =
            tf_motifs::resolve_motif_definition(motif).ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "JASPAR entry '{}' was not found in the local motif registry",
                    motif
                ),
            })?;
        let (llr_matrix, true_log_odds_matrix) =
            Self::prepare_scoring_matrices(&motif_definition.matrix_counts);
        let motif_length_bp = llr_matrix.len();
        if motif_length_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::Internal,
                message: format!(
                    "JASPAR entry '{}' has no usable scoring columns",
                    motif_definition.id
                ),
            });
        }
        if random_sequence_length_bp < motif_length_bp {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "InspectJasparEntry requires random_sequence_length_bp ({}) >= motif length ({}) for '{}'",
                    random_sequence_length_bp, motif_length_bp, motif_definition.id
                ),
            });
        }

        let random_background =
            Self::deterministic_random_dna_bytes(random_sequence_length_bp, random_seed);
        let mut warnings = vec![];
        let remote_metadata = if include_remote_metadata {
            let cached_rows = match Self::load_jaspar_remote_metadata_snapshot_rows(snapshot_path) {
                Ok(rows) => rows,
                Err(err) => {
                    warnings.push(err);
                    BTreeMap::new()
                }
            };
            let cached_metadata = cached_rows
                .get(&motif_definition.id)
                .map(|row| row.remote_metadata.clone());
            if refresh_remote_metadata {
                match self.fetch_jaspar_remote_metadata(&motif_definition.id) {
                    Ok(metadata) => Some(metadata),
                    Err(err) => {
                        warnings.push(err.message);
                        cached_metadata
                    }
                }
            } else {
                cached_metadata
            }
        } else {
            None
        };

        Ok(JasparEntryExpertView {
            schema: JASPAR_ENTRY_EXPERT_VIEW_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            motif_id: motif_definition.id,
            motif_name: motif_definition.name,
            consensus_iupac: motif_definition.consensus_iupac,
            motif_length_bp,
            registry_entry_count: tf_motifs::all_motif_ids().len(),
            requested_token: motif.to_string(),
            random_sequence_length_bp,
            random_seed,
            background_model: "uniform_acgt_lcg".to_string(),
            include_remote_metadata,
            remote_metadata,
            columns: Self::jaspar_expert_columns(&motif_definition.matrix_counts),
            score_panels: vec![
                Self::jaspar_score_distribution_panel(
                    TfbsScoreTrackValueKind::LlrBits,
                    "LLR bits",
                    &llr_matrix,
                    &random_background,
                ),
                Self::jaspar_score_distribution_panel(
                    TfbsScoreTrackValueKind::TrueLogOddsBits,
                    "True log-odds bits",
                    &true_log_odds_matrix,
                    &random_background,
                ),
            ],
            warnings,
        })
    }

    pub(crate) fn list_jaspar_catalog(
        &self,
        filter: Option<&str>,
        limit: Option<usize>,
        include_remote_metadata: bool,
        refresh_remote_metadata: bool,
    ) -> Result<JasparCatalogReport, EngineError> {
        self.list_jaspar_catalog_with_snapshot_path(
            filter,
            limit,
            include_remote_metadata,
            refresh_remote_metadata,
            None,
        )
    }

    pub(crate) fn list_jaspar_catalog_with_snapshot_path(
        &self,
        filter: Option<&str>,
        limit: Option<usize>,
        include_remote_metadata: bool,
        refresh_remote_metadata: bool,
        snapshot_path: Option<&str>,
    ) -> Result<JasparCatalogReport, EngineError> {
        let requested_filter = filter.map(str::trim).filter(|value| !value.is_empty());
        let (registry_entry_count, _, selected_rows) =
            self.select_jaspar_catalog_rows(&[], requested_filter, limit)?;
        let mut warnings = vec![];
        let cached_rows = if include_remote_metadata {
            match Self::load_jaspar_remote_metadata_snapshot_rows(snapshot_path) {
                Ok(rows) => rows,
                Err(err) => {
                    warnings.push(err);
                    BTreeMap::new()
                }
            }
        } else {
            BTreeMap::new()
        };
        let mut rows = vec![];
        for row in selected_rows {
            let remote_summary = if include_remote_metadata {
                let cached_summary = cached_rows.get(&row.id).and_then(|entry| {
                    entry.remote_summary.clone().or_else(|| {
                        Some(Self::jaspar_remote_metadata_summary(&entry.remote_metadata))
                    })
                });
                if refresh_remote_metadata {
                    match self.fetch_jaspar_remote_metadata(&row.id) {
                        Ok(metadata) => Some(Self::jaspar_remote_metadata_summary(&metadata)),
                        Err(err) => {
                            warnings.push(format!("{}: {}", row.id, err.message));
                            cached_summary
                        }
                    }
                } else {
                    cached_summary
                }
            } else {
                None
            };
            rows.push(JasparCatalogRow {
                motif_id: row.id,
                motif_name: row.name,
                consensus_iupac: row.consensus_iupac,
                motif_length_bp: row.motif_length_bp,
                remote_summary,
            });
        }

        Ok(JasparCatalogReport {
            schema: JASPAR_CATALOG_REPORT_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            filter: requested_filter.map(str::to_string),
            limit,
            include_remote_metadata,
            registry_entry_count,
            returned_entry_count: rows.len(),
            rows,
            warnings,
        })
    }

    pub(crate) fn summarize_jaspar_entries(
        &self,
        motifs: &[String],
        random_sequence_length_bp: usize,
        random_seed: u64,
    ) -> Result<JasparEntryPresentationReport, EngineError> {
        if random_sequence_length_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SummarizeJasparEntries requires random_sequence_length_bp >= 1"
                    .to_string(),
            });
        }

        let requested_motifs = if motifs.is_empty() {
            tf_motifs::all_motif_ids()
        } else {
            motifs
                .iter()
                .map(|motif| motif.trim().to_string())
                .filter(|motif| !motif.is_empty())
                .collect::<Vec<_>>()
        };
        if requested_motifs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SummarizeJasparEntries requires at least one motif or a non-empty local JASPAR registry".to_string(),
            });
        }

        let random_background =
            Self::deterministic_random_dna_bytes(random_sequence_length_bp, random_seed);

        let mut seen_ids = BTreeSet::new();
        let mut rows = vec![];
        for token in &requested_motifs {
            let motif = tf_motifs::resolve_motif_definition(token).ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "JASPAR entry '{}' was not found in the local motif registry",
                    token
                ),
            })?;
            if !seen_ids.insert(motif.id.clone()) {
                continue;
            }

            let (llr_matrix, true_log_odds_matrix) =
                Self::prepare_scoring_matrices(&motif.matrix_counts);
            if llr_matrix.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::Internal,
                    message: format!("JASPAR entry '{}' has no usable scoring columns", motif.id),
                });
            }
            if random_background.len() < llr_matrix.len() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "random_sequence_length_bp ({}) must be >= motif length ({}) for '{}'",
                        random_background.len(),
                        llr_matrix.len(),
                        motif.id
                    ),
                });
            }

            let llr_panel = Self::jaspar_score_distribution_panel(
                TfbsScoreTrackValueKind::LlrBits,
                "LLR bits",
                &llr_matrix,
                &random_background,
            );
            let true_log_odds_panel = Self::jaspar_score_distribution_panel(
                TfbsScoreTrackValueKind::TrueLogOddsBits,
                "True log-odds bits",
                &true_log_odds_matrix,
                &random_background,
            );

            rows.push(JasparEntryPresentationRow {
                motif_id: motif.id,
                motif_name: motif.name,
                consensus_iupac: motif.consensus_iupac,
                motif_length_bp: llr_matrix.len(),
                maximizing_sequence: llr_panel.maximizing_sequence.clone(),
                minimizing_sequence: llr_panel.minimizing_sequence.clone(),
                maximizing_llr_bits: llr_panel.maximizing_score,
                maximizing_llr_quantile: llr_panel.maximizing_quantile,
                minimizing_llr_bits: llr_panel.minimizing_score,
                minimizing_llr_quantile: llr_panel.minimizing_quantile,
                maximizing_true_log_odds_bits: true_log_odds_panel.maximizing_score,
                maximizing_true_log_odds_quantile: true_log_odds_panel.maximizing_quantile,
                minimizing_true_log_odds_bits: true_log_odds_panel.minimizing_score,
                minimizing_true_log_odds_quantile: true_log_odds_panel.minimizing_quantile,
                llr_bits_distribution: llr_panel.distribution,
                true_log_odds_bits_distribution: true_log_odds_panel.distribution,
            });
        }

        rows.sort_by(|left, right| {
            left.motif_id
                .to_ascii_uppercase()
                .cmp(&right.motif_id.to_ascii_uppercase())
                .then_with(|| {
                    left.motif_name
                        .as_deref()
                        .unwrap_or("")
                        .to_ascii_uppercase()
                        .cmp(
                            &right
                                .motif_name
                                .as_deref()
                                .unwrap_or("")
                                .to_ascii_uppercase(),
                        )
                })
        });

        Ok(JasparEntryPresentationReport {
            schema: JASPAR_ENTRY_PRESENTATION_REPORT_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            requested_motifs,
            registry_entry_count: tf_motifs::all_motif_ids().len(),
            resolved_entry_count: rows.len(),
            random_sequence_length_bp,
            random_seed,
            background_model: "uniform_acgt_lcg".to_string(),
            rows,
        })
    }

    pub(crate) fn write_jaspar_entry_presentation_report_json(
        &self,
        report: &JasparEntryPresentationReport,
        path: &str,
    ) -> Result<(), EngineError> {
        self.write_pretty_json_file(report, path, "JASPAR entry presentation report")
    }

    pub(crate) fn write_jaspar_catalog_report_json(
        &self,
        report: &JasparCatalogReport,
        path: &str,
    ) -> Result<(), EngineError> {
        self.write_pretty_json_file(report, path, "JASPAR catalog report")
    }

    pub(crate) fn write_jaspar_remote_metadata_snapshot_json(
        &self,
        snapshot: &JasparRemoteMetadataSnapshot,
        path: &str,
    ) -> Result<(), EngineError> {
        self.write_pretty_json_file(snapshot, path, "JASPAR remote-metadata snapshot")
    }

    pub(crate) fn write_jaspar_entry_expert_view_json(
        &self,
        report: &JasparEntryExpertView,
        path: &str,
    ) -> Result<(), EngineError> {
        self.write_pretty_json_file(report, path, "JASPAR entry expert view")
    }
}
