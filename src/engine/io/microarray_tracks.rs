//! Microarray track manifest parsing and genome-anchor projection helpers.
//!
//! This keeps Clariom D array coordinate checks and feature materialization in
//! the shared engine so GUI, CLI, and automation paths all reject mismatched
//! genome builds the same way.

use super::*;

impl GentleEngine {
    /// Load and validate a microarray track manifest without projecting rows.
    pub fn inspect_microarray_track_manifest(
        &self,
        manifest_path: &str,
    ) -> Result<MicroarrayTrackManifest, EngineError> {
        Self::load_microarray_track_manifest(manifest_path)
    }

    pub(super) fn is_generated_microarray_track_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated")
            .any(|v| v.eq_ignore_ascii_case(MICROARRAY_TRACK_GENERATED_TAG))
            || feature
                .qualifier_values("gentle_track_source")
                .any(|v| v.eq_ignore_ascii_case("Array"))
    }

    pub(super) fn remove_generated_microarray_track_features(
        features: &mut Vec<gb_io::seq::Feature>,
    ) {
        features.retain(|feature| !Self::is_generated_microarray_track_feature(feature));
    }

    pub(super) fn load_microarray_track_manifest(
        manifest_path: &str,
    ) -> Result<MicroarrayTrackManifest, EngineError> {
        let path = manifest_path.trim();
        if path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Microarray track manifest path must not be empty".to_string(),

                cause_chain: vec![],
            });
        }
        let file = File::open(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not open microarray track manifest '{path}': {e}"),

            cause_chain: vec![],
        })?;
        let mut manifest: MicroarrayTrackManifest = serde_json::from_reader(BufReader::new(file))
            .map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not parse microarray track manifest '{path}': {e}"),

            cause_chain: vec![],
        })?;
        if manifest.schema.trim() != MICROARRAY_TRACK_MANIFEST_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Microarray track manifest '{}' uses schema '{}' but '{}' is required",
                    path,
                    manifest.schema.trim(),
                    MICROARRAY_TRACK_MANIFEST_SCHEMA
                ),

                cause_chain: vec![],
            });
        }
        manifest.dataset = manifest.dataset.trim().to_string();
        manifest.platform = manifest.platform.trim().to_string();
        manifest.normalization = manifest.normalization.trim().to_string();
        manifest.coordinate_system = manifest.coordinate_system.trim().to_string();
        manifest.supported_genome_ids =
            Self::dedupe_nonempty_strings(&manifest.supported_genome_ids);
        manifest.contrast_order = Self::dedupe_nonempty_strings(&manifest.contrast_order);
        for contrast in &mut manifest.contrasts {
            contrast.contrast = contrast.contrast.trim().to_string();
            contrast.level = contrast.level.trim().to_string();
            contrast.path = contrast.path.trim().to_string();
        }
        manifest.contrasts.retain(|contrast| {
            !contrast.contrast.is_empty() && !contrast.level.is_empty() && !contrast.path.is_empty()
        });
        for projection in &mut manifest.coordinate_projections {
            projection.source_genome_id = projection.source_genome_id.trim().to_string();
            projection.target_genome_id = projection.target_genome_id.trim().to_string();
            projection.method = projection.method.trim().to_string();
            projection.path = projection.path.trim().to_string();
        }
        manifest.coordinate_projections.retain(|projection| {
            !projection.source_genome_id.is_empty()
                && !projection.target_genome_id.is_empty()
                && !projection.path.is_empty()
        });
        if manifest.dataset.is_empty() {
            manifest.dataset = "unknown".to_string();
        }
        if manifest.platform.is_empty() {
            manifest.platform = "unknown".to_string();
        }
        if manifest.normalization.is_empty() {
            manifest.normalization = "unknown".to_string();
        }
        manifest.source_path = Some(path.to_string());
        Ok(manifest)
    }

    fn dedupe_nonempty_strings(values: &[String]) -> Vec<String> {
        let mut out = Vec::new();
        let mut seen = BTreeSet::new();
        for value in values {
            let trimmed = value.trim();
            if trimmed.is_empty() {
                continue;
            }
            let key = trimmed.to_ascii_lowercase();
            if seen.insert(key) {
                out.push(trimmed.to_string());
            }
        }
        out
    }

    fn microarray_manifest_supports_anchor(
        manifest: &MicroarrayTrackManifest,
        anchor: &GenomeSequenceAnchor,
    ) -> bool {
        let anchor_id = anchor.genome_id.trim();
        if anchor_id.is_empty() {
            return false;
        }
        let coordinate = manifest.coordinate_system.trim();
        if !coordinate.is_empty() && coordinate.eq_ignore_ascii_case(anchor_id) {
            return true;
        }
        manifest
            .supported_genome_ids
            .iter()
            .any(|id| id.trim().eq_ignore_ascii_case(anchor_id))
    }

    fn normalize_genome_build_token(value: &str) -> String {
        value
            .chars()
            .filter(|ch| ch.is_ascii_alphanumeric())
            .flat_map(|ch| ch.to_lowercase())
            .collect()
    }

    fn genome_build_tokens_match(left: &str, right: &str) -> bool {
        let left = Self::normalize_genome_build_token(left);
        let right = Self::normalize_genome_build_token(right);
        !left.is_empty() && left == right
    }

    fn resolve_coordinate_projection_spec<'a>(
        manifest: &'a MicroarrayTrackManifest,
        anchor: &GenomeSequenceAnchor,
    ) -> Option<&'a GenomeCoordinateProjectionSpec> {
        manifest.coordinate_projections.iter().find(|projection| {
            Self::genome_build_tokens_match(
                &projection.source_genome_id,
                &manifest.coordinate_system,
            ) && Self::genome_build_tokens_match(&projection.target_genome_id, &anchor.genome_id)
        })
    }

    fn parse_projection_optional_strand(value: Option<&str>) -> Option<char> {
        value
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .and_then(|value| value.chars().next())
            .filter(|strand| matches!(strand, '+' | '-'))
    }

    fn parse_projection_usize(value: Option<&str>, label: &str) -> Result<usize, String> {
        value
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .ok_or_else(|| format!("missing {label}"))?
            .parse::<usize>()
            .map_err(|_| format!("invalid {label}"))
    }

    fn parse_genome_projection_block(
        line: &str,
        header: &HashMap<String, usize>,
        default_method: &str,
    ) -> Result<GenomeCoordinateProjectionBlock, String> {
        let fields = line.split('\t').collect::<Vec<_>>();
        let source_genome_id = Self::microarray_tsv_field(
            &fields,
            header,
            &["source_genome_id", "source_genome", "from_genome_id"],
        )
        .unwrap_or("")
        .trim()
        .to_string();
        let target_genome_id = Self::microarray_tsv_field(
            &fields,
            header,
            &["target_genome_id", "target_genome", "to_genome_id"],
        )
        .unwrap_or("")
        .trim()
        .to_string();
        let source_chrom = Self::microarray_tsv_field(
            &fields,
            header,
            &["source_chrom", "source_chromosome", "chrom", "chromosome"],
        )
        .unwrap_or("")
        .trim()
        .to_string();
        let target_chrom =
            Self::microarray_tsv_field(&fields, header, &["target_chrom", "target_chromosome"])
                .unwrap_or("")
                .trim()
                .to_string();
        if source_genome_id.is_empty() {
            return Err("missing source_genome_id".to_string());
        }
        if target_genome_id.is_empty() {
            return Err("missing target_genome_id".to_string());
        }
        if source_chrom.is_empty() {
            return Err("missing source_chrom".to_string());
        }
        if target_chrom.is_empty() {
            return Err("missing target_chrom".to_string());
        }
        let source_start_1based = Self::parse_projection_usize(
            Self::microarray_tsv_field(&fields, header, &["source_start_1based", "source_start"]),
            "source_start_1based",
        )?;
        let source_end_1based = Self::parse_projection_usize(
            Self::microarray_tsv_field(&fields, header, &["source_end_1based", "source_end"]),
            "source_end_1based",
        )?;
        let target_start_1based = Self::parse_projection_usize(
            Self::microarray_tsv_field(&fields, header, &["target_start_1based", "target_start"]),
            "target_start_1based",
        )?;
        let target_end_1based = Self::parse_projection_usize(
            Self::microarray_tsv_field(&fields, header, &["target_end_1based", "target_end"]),
            "target_end_1based",
        )?;
        if source_start_1based == 0
            || target_start_1based == 0
            || source_end_1based < source_start_1based
            || target_end_1based < target_start_1based
        {
            return Err("invalid projection interval".to_string());
        }
        let method = Self::microarray_tsv_field(&fields, header, &["method"])
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or(default_method)
            .to_string();
        Ok(GenomeCoordinateProjectionBlock {
            source_genome_id,
            target_genome_id,
            source_chrom,
            source_start_1based,
            source_end_1based,
            source_strand: Self::parse_projection_optional_strand(Self::microarray_tsv_field(
                &fields,
                header,
                &["source_strand"],
            )),
            target_chrom,
            target_start_1based,
            target_end_1based,
            target_strand: Self::parse_projection_optional_strand(Self::microarray_tsv_field(
                &fields,
                header,
                &["target_strand"],
            )),
            method,
        })
    }

    fn load_genome_coordinate_projection_blocks(
        projection_path: &str,
        default_method: &str,
        source_genome_id: &str,
        target_genome_id: &str,
    ) -> Result<Vec<GenomeCoordinateProjectionBlock>, EngineError> {
        let mut reader = Self::open_text_reader(projection_path)?;
        let mut line = String::new();
        let mut header: Option<HashMap<String, usize>> = None;
        let mut blocks = Vec::new();
        let mut line_no = 0usize;
        while {
            line.clear();
            reader.read_line(&mut line).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not read genome projection map '{projection_path}': {e}"),

                cause_chain: vec![],
            })? > 0
        } {
            line_no += 1;
            let trimmed = line.trim_end_matches(['\r', '\n']);
            if trimmed.trim().is_empty() || trimmed.trim_start().starts_with('#') {
                continue;
            }
            if header.is_none() {
                header = Some(
                    trimmed
                        .split('\t')
                        .enumerate()
                        .map(|(idx, key)| (Self::normalized_microarray_header_key(key), idx))
                        .collect(),
                );
                continue;
            }
            let block = Self::parse_genome_projection_block(
                trimmed,
                header.as_ref().expect("header set"),
                default_method,
            )
            .map_err(|err| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Genome projection map '{}' line {} skipped: {}",
                    projection_path, line_no, err
                ),

                cause_chain: vec![],
            })?;
            if Self::genome_build_tokens_match(&block.source_genome_id, source_genome_id)
                && Self::genome_build_tokens_match(&block.target_genome_id, target_genome_id)
            {
                blocks.push(block);
            }
        }
        if blocks.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Genome projection map '{}' contains no blocks for {} -> {}",
                    projection_path, source_genome_id, target_genome_id
                ),

                cause_chain: vec![],
            });
        }
        Ok(blocks)
    }

    fn project_genome_interval_with_blocks(
        blocks: &[GenomeCoordinateProjectionBlock],
        chrom: &str,
        start_1based: usize,
        end_1based: usize,
        strand: Option<char>,
    ) -> Option<ProjectedGenomeInterval> {
        let block = blocks.iter().find(|block| {
            Self::chromosomes_match(&block.source_chrom, chrom)
                && start_1based >= block.source_start_1based
                && end_1based <= block.source_end_1based
        })?;
        let source_len = block.source_end_1based - block.source_start_1based + 1;
        let target_len = block.target_end_1based - block.target_start_1based + 1;
        if source_len != target_len {
            return None;
        }
        let left_offset = start_1based - block.source_start_1based;
        let right_offset = end_1based - block.source_start_1based;
        let target_reverse = block.target_strand == Some('-');
        let (target_start_1based, target_end_1based) = if target_reverse {
            (
                block.target_end_1based.saturating_sub(right_offset),
                block.target_end_1based.saturating_sub(left_offset),
            )
        } else {
            (
                block.target_start_1based + left_offset,
                block.target_start_1based + right_offset,
            )
        };
        let target_strand = match (strand, target_reverse) {
            (Some('+'), true) => Some('-'),
            (Some('-'), true) => Some('+'),
            (Some(value), false) => Some(value),
            _ => block.source_strand.or(strand),
        };
        Some(ProjectedGenomeInterval {
            target_chrom: block.target_chrom.clone(),
            target_start_1based,
            target_end_1based,
            target_strand,
            method: block.method.clone(),
            status: "mapped_unique_interval_block".to_string(),
        })
    }

    pub(super) fn project_genome_interval_from_map(
        source_genome_id: &str,
        target_genome_id: &str,
        projection_path: &str,
        method: &str,
        chrom: &str,
        start_1based: usize,
        end_1based: usize,
        strand: Option<char>,
    ) -> Result<GenomeCoordinateProjectionReport, EngineError> {
        if start_1based == 0 || end_1based < start_1based {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ProjectGenomeInterval requires a valid 1-based inclusive interval"
                    .to_string(),

                cause_chain: vec![],
            });
        }
        let blocks = Self::load_genome_coordinate_projection_blocks(
            projection_path,
            method,
            source_genome_id,
            target_genome_id,
        )?;
        let projected = Self::project_genome_interval_with_blocks(
            &blocks,
            chrom,
            start_1based,
            end_1based,
            strand,
        );
        let mut report = GenomeCoordinateProjectionReport {
            schema: GENOME_COORDINATE_PROJECTION_REPORT_SCHEMA.to_string(),
            source_genome_id: source_genome_id.to_string(),
            target_genome_id: target_genome_id.to_string(),
            projection_path: projection_path.to_string(),
            method: method.to_string(),
            input_chrom: chrom.to_string(),
            input_start_1based: start_1based,
            input_end_1based: end_1based,
            input_strand: strand.map(|value| value.to_string()),
            mapped: projected.is_some(),
            ..Default::default()
        };
        if let Some(projected) = projected {
            report.intervals.push(GenomeCoordinateProjectionInterval {
                source_chrom: chrom.to_string(),
                source_start_1based: start_1based,
                source_end_1based: end_1based,
                source_strand: strand.map(|value| value.to_string()),
                target_chrom: projected.target_chrom,
                target_start_1based: projected.target_start_1based,
                target_end_1based: projected.target_end_1based,
                target_strand: projected.target_strand.map(|value| value.to_string()),
                method: projected.method,
                status: projected.status,
            });
        } else {
            report.warnings.push(format!(
                "No projection block mapped {}:{}-{} from {} to {}",
                chrom, start_1based, end_1based, source_genome_id, target_genome_id
            ));
        }
        Ok(report)
    }

    fn resolve_microarray_contrasts(
        manifest: &MicroarrayTrackManifest,
        requested_contrasts: &[String],
        level: &str,
    ) -> Result<Vec<MicroarrayTrackContrast>, EngineError> {
        let level = level.trim();
        if level.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ProjectMicroarrayTrack level must not be empty".to_string(),

                cause_chain: vec![],
            });
        }
        let available = manifest
            .contrasts
            .iter()
            .filter(|entry| entry.level.eq_ignore_ascii_case(level))
            .cloned()
            .collect::<Vec<_>>();
        if available.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Microarray track manifest has no '{}' contrast TSVs", level),

                cause_chain: vec![],
            });
        }

        let requested = requested_contrasts
            .iter()
            .map(|contrast| contrast.trim())
            .filter(|contrast| !contrast.is_empty())
            .map(|contrast| contrast.to_string())
            .collect::<Vec<_>>();
        let ordered_names = if requested.is_empty() {
            if manifest.contrast_order.is_empty() {
                available
                    .iter()
                    .map(|entry| entry.contrast.clone())
                    .collect::<Vec<_>>()
            } else {
                manifest.contrast_order.clone()
            }
        } else {
            requested
        };

        let mut out = Vec::new();
        let mut missing = Vec::new();
        for name in ordered_names {
            if let Some(entry) = available
                .iter()
                .find(|entry| entry.contrast.eq_ignore_ascii_case(&name))
            {
                if !out.iter().any(|existing: &MicroarrayTrackContrast| {
                    existing.contrast.eq_ignore_ascii_case(&entry.contrast)
                }) {
                    out.push(entry.clone());
                }
            } else {
                missing.push(name);
            }
        }
        if !missing.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Microarray track manifest lacks requested contrast(s) for level '{}': {}",
                    level,
                    missing.join(", ")
                ),

                cause_chain: vec![],
            });
        }
        Ok(out)
    }

    fn resolve_microarray_track_path(manifest_path: &str, track_path: &str) -> String {
        let path = Path::new(track_path);
        if path.is_absolute() {
            return path.to_string_lossy().to_string();
        }
        Path::new(manifest_path)
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(path)
            .to_string_lossy()
            .to_string()
    }

    fn normalized_microarray_header_key(key: &str) -> String {
        key.chars()
            .filter(|ch| ch.is_ascii_alphanumeric())
            .map(|ch| ch.to_ascii_lowercase())
            .collect()
    }

    fn microarray_tsv_field<'a>(
        fields: &'a [&'a str],
        header: &HashMap<String, usize>,
        aliases: &[&str],
    ) -> Option<&'a str> {
        for alias in aliases {
            let key = Self::normalized_microarray_header_key(alias);
            if let Some(index) = header.get(&key)
                && let Some(value) = fields.get(*index)
            {
                let trimmed = value.trim();
                if !trimmed.is_empty() && trimmed != "." && !trimmed.eq_ignore_ascii_case("NA") {
                    return Some(trimmed);
                }
            }
        }
        None
    }

    fn parse_microarray_optional_f64(value: Option<&str>) -> Option<f64> {
        value.and_then(|raw| {
            let trimmed = raw.trim();
            if trimmed.is_empty()
                || trimmed.eq_ignore_ascii_case("NA")
                || trimmed.eq_ignore_ascii_case("NaN")
                || trimmed == "."
            {
                None
            } else {
                trimmed
                    .parse::<f64>()
                    .ok()
                    .filter(|value| value.is_finite())
            }
        })
    }

    fn parse_microarray_track_row(
        line: &str,
        header: &HashMap<String, usize>,
    ) -> Result<MicroarrayTrackRow, String> {
        let fields = line.split('\t').collect::<Vec<_>>();
        let chromosome = Self::microarray_tsv_field(
            &fields,
            header,
            &["chrom", "chromosome", "seqnames", "genoName"],
        )
        .ok_or_else(|| "missing chrom/chromosome".to_string())?
        .to_string();
        let start_raw = Self::microarray_tsv_field(
            &fields,
            header,
            &["start_1based", "genomic_start_1based", "start", "pos_start"],
        )
        .ok_or_else(|| "missing start_1based".to_string())?;
        let end_raw = Self::microarray_tsv_field(
            &fields,
            header,
            &["end_1based", "genomic_end_1based", "end", "stop", "pos_end"],
        )
        .ok_or_else(|| "missing end_1based".to_string())?;
        let start_1based = start_raw
            .parse::<usize>()
            .map_err(|e| format!("invalid start_1based '{start_raw}': {e}"))?;
        let end_1based = end_raw
            .parse::<usize>()
            .map_err(|e| format!("invalid end_1based '{end_raw}': {e}"))?;
        if start_1based == 0 || end_1based < start_1based {
            return Err(format!(
                "invalid microarray interval {chromosome}:{start_1based}-{end_1based}"
            ));
        }
        let strand = Self::microarray_tsv_field(&fields, header, &["strand"])
            .and_then(|value| value.chars().next())
            .filter(|strand| matches!(strand, '+' | '-'));
        let feature_id = Self::microarray_tsv_field(
            &fields,
            header,
            &["feature_id", "probeset_id", "probe_set_id", "id"],
        )
        .unwrap_or("")
        .to_string();
        Ok(MicroarrayTrackRow {
            chromosome,
            start_1based,
            end_1based,
            strand,
            feature_id,
            transcript_cluster_id: Self::microarray_tsv_field(
                &fields,
                header,
                &["transcript_cluster_id", "man_fsetid", "transcriptclusterid"],
            )
            .map(str::to_string),
            exon_id: Self::microarray_tsv_field(&fields, header, &["exon_id", "exonid"])
                .map(str::to_string),
            probe_type: Self::microarray_tsv_field(
                &fields,
                header,
                &["probe_type", "probeset_type", "type", "crosshyb_type"],
            )
            .map(str::to_string),
            logfc: Self::parse_microarray_optional_f64(Self::microarray_tsv_field(
                &fields,
                header,
                &["logFC", "log_fc", "log2FC"],
            )),
            ave_expr: Self::parse_microarray_optional_f64(Self::microarray_tsv_field(
                &fields,
                header,
                &["AveExpr", "average_expression", "ave_expr"],
            )),
            p_value: Self::parse_microarray_optional_f64(Self::microarray_tsv_field(
                &fields,
                header,
                &["P.Value", "P_Value", "pvalue", "p_value"],
            )),
            adj_p_value: Self::parse_microarray_optional_f64(Self::microarray_tsv_field(
                &fields,
                header,
                &["adj.P.Val", "adj_P_Val", "adj_p_value", "FDR"],
            )),
            junction_start_edge: Self::microarray_tsv_field(
                &fields,
                header,
                &["junction_start_edge"],
            )
            .map(str::to_string),
            junction_stop_edge: Self::microarray_tsv_field(
                &fields,
                header,
                &["junction_stop_edge"],
            )
            .map(str::to_string),
            junction_sequence: Self::microarray_tsv_field(&fields, header, &["junction_sequence"])
                .map(str::to_string),
            gene_symbol: Self::microarray_tsv_field(
                &fields,
                header,
                &["gene_symbol", "gene", "gene_assignment"],
            )
            .map(str::to_string),
            has_cds: Self::microarray_tsv_field(&fields, header, &["has_cds"]).map(str::to_string),
        })
    }

    fn format_microarray_number(value: f64) -> String {
        let abs = value.abs();
        if abs > 0.0 && !(0.0001..100_000.0).contains(&abs) {
            format!("{value:.6e}")
        } else {
            format!("{value:.6}")
        }
    }

    fn microarray_feature_group_key(row: &MicroarrayTrackRow) -> String {
        if !row.feature_id.trim().is_empty() {
            return row.feature_id.trim().to_string();
        }
        let transcript = row.transcript_cluster_id.as_deref().unwrap_or("-");
        let exon = row.exon_id.as_deref().unwrap_or("-");
        format!(
            "{}:{}-{}:{transcript}:{exon}",
            row.chromosome, row.start_1based, row.end_1based
        )
    }

    fn build_microarray_feature(
        manifest: &MicroarrayTrackManifest,
        anchor: &GenomeSequenceAnchor,
        manifest_path: &str,
        contrast: &str,
        level: &str,
        row: &MicroarrayTrackRow,
        native_row: Option<&MicroarrayTrackRow>,
        projection: Option<&ProjectedGenomeInterval>,
        local_start_0based: usize,
        local_end_0based_exclusive: usize,
        local_strand: Option<char>,
    ) -> gb_io::seq::Feature {
        let native_row = native_row.unwrap_or(row);
        let group_key = Self::microarray_feature_group_key(native_row);
        let label = if group_key.trim().is_empty() {
            contrast.to_string()
        } else {
            format!("{contrast} {group_key}")
        };
        let assembly_check = if projection.is_some() {
            "projected_from_native_coordinate_system"
        } else if manifest
            .coordinate_system
            .trim()
            .eq_ignore_ascii_case(anchor.genome_id.trim())
        {
            "coordinate_system_matches_anchor"
        } else {
            "supported_genome_id_alias_matches_anchor"
        };
        let mut qualifiers = vec![
            ("label".into(), Some(label)),
            (
                "note".into(),
                Some(format!(
                    "Array track '{}' from '{}' [{}:{}-{}]",
                    contrast, manifest_path, row.chromosome, row.start_1based, row.end_1based
                )),
            ),
            ("gentle_track_source".into(), Some("Array".to_string())),
            (
                "gentle_generated".into(),
                Some(MICROARRAY_TRACK_GENERATED_TAG.to_string()),
            ),
            (
                "gentle_track_name".into(),
                Some(format!("{} {}", manifest.platform, contrast)),
            ),
            ("gentle_track_file".into(), Some(manifest_path.to_string())),
            (
                "gentle_array_dataset".into(),
                Some(manifest.dataset.clone()),
            ),
            (
                "gentle_array_platform".into(),
                Some(manifest.platform.clone()),
            ),
            (
                "gentle_array_normalization".into(),
                Some(manifest.normalization.clone()),
            ),
            (
                "gentle_array_coordinate_system".into(),
                Some(manifest.coordinate_system.clone()),
            ),
            (
                "gentle_array_supported_genome_ids".into(),
                Some(manifest.supported_genome_ids.join(",")),
            ),
            (
                "gentle_array_anchor_genome_id".into(),
                Some(anchor.genome_id.clone()),
            ),
            (
                "gentle_array_anchor_chromosome".into(),
                Some(anchor.chromosome.clone()),
            ),
            (
                "gentle_array_anchor_start_1based".into(),
                Some(anchor.start_1based.to_string()),
            ),
            (
                "gentle_array_anchor_end_1based".into(),
                Some(anchor.end_1based.to_string()),
            ),
            (
                "gentle_array_anchor_strand".into(),
                Some(anchor.strand.unwrap_or('+').to_string()),
            ),
            (
                "gentle_array_assembly_check".into(),
                Some(assembly_check.to_string()),
            ),
            ("gentle_array_contrast".into(), Some(contrast.to_string())),
            ("gentle_array_level".into(), Some(level.to_string())),
            ("chromosome".into(), Some(row.chromosome.clone())),
            ("start_1based".into(), Some(row.start_1based.to_string())),
            ("end_1based".into(), Some(row.end_1based.to_string())),
            (
                "genomic_start_1based".into(),
                Some(row.start_1based.to_string()),
            ),
            (
                "genomic_end_1based".into(),
                Some(row.end_1based.to_string()),
            ),
            (
                "gentle_array_native_chromosome".into(),
                Some(native_row.chromosome.clone()),
            ),
            (
                "gentle_array_native_start_1based".into(),
                Some(native_row.start_1based.to_string()),
            ),
            (
                "gentle_array_native_end_1based".into(),
                Some(native_row.end_1based.to_string()),
            ),
        ];
        if let Some(projection) = projection {
            qualifiers.push((
                "gentle_array_projection_method".into(),
                Some(projection.method.clone()),
            ));
            qualifiers.push((
                "gentle_array_projection_status".into(),
                Some(projection.status.clone()),
            ));
        }
        if let Some(strand) = native_row.strand {
            qualifiers.push(("array_strand".into(), Some(strand.to_string())));
        }
        if let Some(strand) = local_strand {
            qualifiers.push(("strand".into(), Some(strand.to_string())));
        }
        if !row.feature_id.trim().is_empty() {
            qualifiers.push(("feature_id".into(), Some(row.feature_id.clone())));
            qualifiers.push((
                "gentle_array_feature_id".into(),
                Some(row.feature_id.clone()),
            ));
        }
        if let Some(value) = &row.transcript_cluster_id {
            qualifiers.push(("transcript_cluster_id".into(), Some(value.clone())));
        }
        if let Some(value) = &row.exon_id {
            qualifiers.push(("exon_id".into(), Some(value.clone())));
        }
        if let Some(value) = &row.probe_type {
            qualifiers.push(("gentle_array_probe_type".into(), Some(value.clone())));
        }
        if let Some(value) = row.logfc {
            let formatted = Self::format_microarray_number(value);
            qualifiers.push(("logFC".into(), Some(formatted.clone())));
            qualifiers.push(("score".into(), Some(formatted)));
        }
        if let Some(value) = row.ave_expr {
            qualifiers.push((
                "AveExpr".into(),
                Some(Self::format_microarray_number(value)),
            ));
        }
        if let Some(value) = row.p_value {
            qualifiers.push((
                "P_Value".into(),
                Some(Self::format_microarray_number(value)),
            ));
        }
        if let Some(value) = row.adj_p_value {
            qualifiers.push((
                "adj_P_Val".into(),
                Some(Self::format_microarray_number(value)),
            ));
        }
        if let Some(value) = &row.junction_start_edge {
            qualifiers.push(("junction_start_edge".into(), Some(value.clone())));
        }
        if let Some(value) = &row.junction_stop_edge {
            qualifiers.push(("junction_stop_edge".into(), Some(value.clone())));
        }
        if let Some(value) = &row.junction_sequence {
            qualifiers.push(("junction_sequence".into(), Some(value.clone())));
        }
        if let Some(value) = &row.gene_symbol {
            qualifiers.push(("gene".into(), Some(value.clone())));
            qualifiers.push(("gene_symbol".into(), Some(value.clone())));
        }
        if let Some(value) = &row.has_cds {
            qualifiers.push(("has_cds".into(), Some(value.clone())));
        }

        let base_location = gb_io::seq::Location::simple_range(
            local_start_0based as i64,
            local_end_0based_exclusive as i64,
        );
        let location = if local_strand == Some('-') {
            gb_io::seq::Location::Complement(Box::new(base_location))
        } else {
            base_location
        };

        gb_io::seq::Feature {
            kind: "track".into(),
            location,
            qualifiers,
        }
    }

    pub(super) fn project_microarray_track(
        dna: &mut DNAsequence,
        anchor: &GenomeSequenceAnchor,
        seq_id: &str,
        manifest_path: &str,
        contrasts: &[String],
        level: &str,
        min_abs_logfc: Option<f64>,
        max_adj_p: Option<f64>,
        max_features: Option<usize>,
        clear_existing: bool,
    ) -> Result<MicroarrayProjectionReport, EngineError> {
        let manifest = Self::load_microarray_track_manifest(manifest_path)?;
        let direct_coordinate_match = Self::microarray_manifest_supports_anchor(&manifest, anchor);
        let projection_spec = if direct_coordinate_match {
            None
        } else {
            Self::resolve_coordinate_projection_spec(&manifest, anchor)
        };
        let projection_path = projection_spec
            .map(|spec| Self::resolve_microarray_track_path(manifest_path, &spec.path));
        let projection_blocks =
            if let (Some(spec), Some(path)) = (projection_spec, &projection_path) {
                Some(Self::load_genome_coordinate_projection_blocks(
                    path,
                    &spec.method,
                    &spec.source_genome_id,
                    &spec.target_genome_id,
                )?)
            } else {
                None
            };
        if !direct_coordinate_match && projection_blocks.is_none() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Microarray track manifest '{}' coordinate system '{}' is not declared compatible with sequence anchor genome_id '{}' (supported: {}; projection maps: {})",
                    manifest_path,
                    manifest.coordinate_system,
                    anchor.genome_id,
                    if manifest.supported_genome_ids.is_empty() {
                        "none".to_string()
                    } else {
                        manifest.supported_genome_ids.join(", ")
                    },
                    if manifest.coordinate_projections.is_empty() {
                        "none".to_string()
                    } else {
                        manifest
                            .coordinate_projections
                            .iter()
                            .map(|projection| {
                                format!(
                                    "{}->{}",
                                    projection.source_genome_id, projection.target_genome_id
                                )
                            })
                            .collect::<Vec<_>>()
                            .join(", ")
                    }
                ),

                cause_chain: vec![],
            });
        }
        let level = if level.trim().is_empty() {
            "probeset"
        } else {
            level.trim()
        };
        let selected_contrasts = Self::resolve_microarray_contrasts(&manifest, contrasts, level)?;
        let feature_limit = max_features.unwrap_or(MAX_IMPORTED_SIGNAL_FEATURES);
        if feature_limit == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ProjectMicroarrayTrack max_features must be greater than zero"
                    .to_string(),

                cause_chain: vec![],
            });
        }
        if let Some(min_abs_logfc) = min_abs_logfc
            && (!min_abs_logfc.is_finite() || min_abs_logfc < 0.0)
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ProjectMicroarrayTrack min_abs_logfc must be >= 0".to_string(),

                cause_chain: vec![],
            });
        }
        if let Some(max_adj_p) = max_adj_p
            && (!max_adj_p.is_finite() || !(0.0..=1.0).contains(&max_adj_p))
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ProjectMicroarrayTrack max_adj_p must be within 0..=1".to_string(),

                cause_chain: vec![],
            });
        }

        let mut report = MicroarrayProjectionReport {
            schema: MICROARRAY_PROJECTION_REPORT_SCHEMA.to_string(),
            seq_id: seq_id.to_string(),
            manifest_path: manifest_path.to_string(),
            dataset: manifest.dataset.clone(),
            platform: manifest.platform.clone(),
            normalization: manifest.normalization.clone(),
            coordinate_system: manifest.coordinate_system.clone(),
            coordinate_projection_used: projection_blocks.is_some(),
            coordinate_projection_method: projection_spec.map(|spec| spec.method.clone()),
            coordinate_projection_path: projection_path.clone(),
            anchor_genome_id: anchor.genome_id.clone(),
            anchor_chromosome: anchor.chromosome.clone(),
            anchor_start_1based: anchor.start_1based,
            anchor_end_1based: anchor.end_1based,
            anchor_strand: anchor.strand.unwrap_or('+').to_string(),
            requested_contrasts: contrasts
                .iter()
                .map(|contrast| contrast.trim())
                .filter(|contrast| !contrast.is_empty())
                .map(str::to_string)
                .collect(),
            projected_contrasts: selected_contrasts
                .iter()
                .map(|entry| entry.contrast.clone())
                .collect(),
            level: level.to_string(),
            warnings: manifest.warnings.clone(),
            ..Default::default()
        };

        let mut pending_features: Vec<PendingMicroarrayFeature> = Vec::new();
        let mut mismatch_counts: HashMap<String, usize> = HashMap::new();

        'contrasts: for contrast_entry in selected_contrasts {
            let tsv_path = Self::resolve_microarray_track_path(manifest_path, &contrast_entry.path);
            let mut reader = Self::open_text_reader(&tsv_path)?;
            let mut line = String::new();
            let mut line_no = 0usize;
            let mut header: Option<HashMap<String, usize>> = None;
            while {
                line.clear();
                reader.read_line(&mut line).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not read microarray track TSV '{tsv_path}': {e}"),

                    cause_chain: vec![],
                })? > 0
            } {
                line_no += 1;
                let trimmed = line.trim_end_matches(['\r', '\n']);
                if trimmed.trim().is_empty() || trimmed.trim_start().starts_with('#') {
                    continue;
                }
                if header.is_none() {
                    let parsed_header = trimmed
                        .split('\t')
                        .enumerate()
                        .map(|(idx, key)| (Self::normalized_microarray_header_key(key), idx))
                        .filter(|(key, _)| !key.is_empty())
                        .collect::<HashMap<_, _>>();
                    header = Some(parsed_header);
                    continue;
                }
                report.parsed_rows += 1;
                let mut row = match Self::parse_microarray_track_row(
                    trimmed,
                    header.as_ref().expect("header set"),
                ) {
                    Ok(row) => row,
                    Err(err) => {
                        report.skipped_rows += 1;
                        report.skipped_invalid += 1;
                        if report.warnings.len() < 24 {
                            report.warnings.push(format!(
                                "Microarray track '{}' line {} skipped: {}",
                                contrast_entry.contrast, line_no, err
                            ));
                        }
                        continue;
                    }
                };
                let native_row = row.clone();

                if min_abs_logfc
                    .map(|threshold| {
                        row.logfc
                            .map(|value| value.abs() < threshold)
                            .unwrap_or(true)
                    })
                    .unwrap_or(false)
                {
                    report.skipped_rows += 1;
                    report.skipped_filter += 1;
                    continue;
                }
                if max_adj_p
                    .map(|threshold| {
                        row.adj_p_value
                            .map(|value| value > threshold)
                            .unwrap_or(true)
                    })
                    .unwrap_or(false)
                {
                    report.skipped_rows += 1;
                    report.skipped_filter += 1;
                    continue;
                }

                let projected_interval = if let Some(blocks) = projection_blocks.as_deref() {
                    match Self::project_genome_interval_with_blocks(
                        blocks,
                        &native_row.chromosome,
                        native_row.start_1based,
                        native_row.end_1based,
                        native_row.strand,
                    ) {
                        Some(projected) => {
                            row.chromosome = projected.target_chrom.clone();
                            row.start_1based = projected.target_start_1based;
                            row.end_1based = projected.target_end_1based;
                            row.strand = projected.target_strand;
                            Some(projected)
                        }
                        None => {
                            report.skipped_rows += 1;
                            report.skipped_projection_unmapped += 1;
                            continue;
                        }
                    }
                } else {
                    None
                };

                if !Self::chromosomes_match(&row.chromosome, &anchor.chromosome) {
                    report.skipped_rows += 1;
                    report.skipped_wrong_chromosome += 1;
                    *mismatch_counts.entry(row.chromosome.clone()).or_insert(0) += 1;
                    continue;
                }
                if row.end_1based < anchor.start_1based || row.start_1based > anchor.end_1based {
                    report.skipped_rows += 1;
                    report.skipped_non_overlap += 1;
                    continue;
                }
                let overlap_start_1based = row.start_1based.max(anchor.start_1based);
                let overlap_end_1based = row.end_1based.min(anchor.end_1based);
                if overlap_end_1based < overlap_start_1based {
                    report.skipped_rows += 1;
                    report.skipped_non_overlap += 1;
                    continue;
                }
                if pending_features.len() >= feature_limit {
                    report.truncated_at_limit = true;
                    break 'contrasts;
                }

                let (local_start_0based, local_end_0based_exclusive) = if anchor.strand == Some('-')
                {
                    (
                        anchor.end_1based.saturating_sub(overlap_end_1based),
                        anchor.end_1based.saturating_sub(overlap_start_1based) + 1,
                    )
                } else {
                    (
                        overlap_start_1based - anchor.start_1based,
                        overlap_end_1based - anchor.start_1based + 1,
                    )
                };
                let local_strand = match (row.strand, anchor.strand) {
                    (Some('+'), Some('-')) => Some('-'),
                    (Some('-'), Some('-')) => Some('+'),
                    (Some(strand), _) => Some(strand),
                    _ => None,
                };
                let feature = Self::build_microarray_feature(
                    &manifest,
                    anchor,
                    manifest_path,
                    &contrast_entry.contrast,
                    level,
                    &row,
                    Some(&native_row),
                    projected_interval.as_ref(),
                    local_start_0based,
                    local_end_0based_exclusive,
                    local_strand,
                );
                pending_features.push(PendingMicroarrayFeature {
                    feature,
                    group_key: Self::microarray_feature_group_key(&native_row),
                    contrast: contrast_entry.contrast.clone(),
                    logfc: row.logfc,
                    adj_p_value: row.adj_p_value,
                });
            }
            if header.is_none() {
                report.warnings.push(format!(
                    "Microarray track TSV '{}' did not contain a header row",
                    tsv_path
                ));
            }
        }

        if !mismatch_counts.is_empty() {
            let mut sorted = mismatch_counts
                .iter()
                .map(|(chromosome, count)| (chromosome.clone(), *count))
                .collect::<Vec<_>>();
            sorted.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
            let seen = sorted
                .iter()
                .take(3)
                .map(|(chromosome, count)| format!("{chromosome} ({count})"))
                .collect::<Vec<_>>()
                .join(", ");
            report.warnings.push(format!(
                "{} microarray row(s) did not match anchor chromosome '{}' (examples: {})",
                report.skipped_wrong_chromosome, anchor.chromosome, seen
            ));
        }

        let mut grouped_values: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for pending in &pending_features {
            let mut parts = Vec::new();
            if let Some(logfc) = pending.logfc {
                parts.push(format!("logFC={}", Self::format_microarray_number(logfc)));
            }
            if let Some(adj_p_value) = pending.adj_p_value {
                parts.push(format!(
                    "adj.P.Val={}",
                    Self::format_microarray_number(adj_p_value)
                ));
            }
            let summary = if parts.is_empty() {
                pending.contrast.clone()
            } else {
                format!("{} {}", pending.contrast, parts.join(" "))
            };
            grouped_values
                .entry(pending.group_key.clone())
                .or_default()
                .push(summary);
        }

        if clear_existing {
            Self::remove_generated_microarray_track_features(dna.features_mut());
        }

        report.imported_features = pending_features.len();
        for mut pending in pending_features {
            if let Some(values) = grouped_values.get(&pending.group_key) {
                pending
                    .feature
                    .qualifiers
                    .push(("gentle_array_value_summary".into(), Some(values.join("; "))));
            }
            dna.features_mut().push(pending.feature);
        }

        if report.truncated_at_limit {
            report.warnings.push(format!(
                "Microarray projection was truncated after {} features (limit={})",
                report.imported_features, feature_limit
            ));
        }

        Ok(report)
    }
}
