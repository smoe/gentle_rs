//! ENCODE SCREEN cCRE overlap projection and feature materialization.
//!
//! Biological interpretation and assembly checks live here so GUI, CLI, MCP,
//! JavaScript, Lua, and future adapters consume the same evidence contract.

use super::*;

const DEFAULT_ENCODE_CCRE_QUERY_LIMIT: usize = 10_000;
const MAX_ENCODE_CCRE_QUERY_LIMIT: usize = 1_000_000;
const ENCODE_CCRE_EVIDENCE_STATEMENT: &str = "Overlap with an ENCODE SCREEN Registry V4 candidate cis-regulatory element is biochemical candidate-regulatory evidence.";

impl GentleEngine {
    fn normalized_encode_ccre_classes(
        requested: &[String],
        available: &[String],
    ) -> Result<Vec<String>, EngineError> {
        let available_by_lower = available
            .iter()
            .map(|value| (value.to_ascii_lowercase(), value.clone()))
            .collect::<BTreeMap<_, _>>();
        let mut normalized = BTreeSet::new();
        for raw in requested {
            let token = raw.trim();
            if token.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "ENCODE cCRE class filters must not be blank".to_string(),
                    cause_chain: vec![],
                });
            }
            let Some(canonical) = available_by_lower.get(&token.to_ascii_lowercase()) else {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "ENCODE cCRE class '{}' is not present in source; available classes: {}",
                        token,
                        available.join(", ")
                    ),
                    cause_chain: vec![],
                });
            };
            normalized.insert(canonical.clone());
        }
        Ok(normalized.into_iter().collect())
    }

    fn encode_ccre_genomic_query_span(
        anchor: &SequenceGenomeAnchorSummary,
        local_start: usize,
        local_end: usize,
    ) -> Result<(u64, u64), EngineError> {
        let anchor_start = anchor.start_1based.saturating_sub(1) as u64;
        let anchor_end = anchor.end_1based as u64;
        let span = if anchor.strand == Some('-') {
            (
                anchor_end.saturating_sub(local_end as u64),
                anchor_end.saturating_sub(local_start as u64),
            )
        } else {
            (
                anchor_start.saturating_add(local_start as u64),
                anchor_start.saturating_add(local_end as u64),
            )
        };
        if span.1 <= span.0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ENCODE cCRE query resolves to an empty genomic interval".to_string(),
                cause_chain: vec![],
            });
        }
        Ok(span)
    }

    fn project_encode_ccre_interval(
        interval: EncodeCcreInterval,
        anchor: &SequenceGenomeAnchorSummary,
        query_genomic_start: u64,
        query_genomic_end: u64,
    ) -> Option<EncodeCcreOverlapRow> {
        let overlap_start = interval.start_0based.max(query_genomic_start);
        let overlap_end = interval.end_0based_exclusive.min(query_genomic_end);
        if overlap_end <= overlap_start {
            return None;
        }
        let anchor_start = anchor.start_1based.saturating_sub(1) as u64;
        let anchor_end = anchor.end_1based as u64;
        let (local_start, local_end) = if anchor.strand == Some('-') {
            (
                anchor_end.saturating_sub(overlap_end),
                anchor_end.saturating_sub(overlap_start),
            )
        } else {
            (
                overlap_start.saturating_sub(anchor_start),
                overlap_end.saturating_sub(anchor_start),
            )
        };
        Some(EncodeCcreOverlapRow {
            clipped: overlap_start != interval.start_0based
                || overlap_end != interval.end_0based_exclusive,
            interval,
            local_start_0based: usize::try_from(local_start).ok()?,
            local_end_0based_exclusive: usize::try_from(local_end).ok()?,
            genomic_start_0based: overlap_start,
            genomic_end_0based_exclusive: overlap_end,
            overlap_bp: usize::try_from(overlap_end.saturating_sub(overlap_start)).ok()?,
        })
    }

    pub fn query_encode_ccre_overlaps(
        &self,
        seq_id: &str,
        index_path: &str,
        bed_path_override: Option<&str>,
        start_0based: Option<usize>,
        end_0based_exclusive: Option<usize>,
        classes: &[String],
        limit: Option<usize>,
    ) -> Result<EncodeCcreOverlapReport, EngineError> {
        let seq_id = seq_id.trim();
        if seq_id.is_empty() || index_path.trim().is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "QueryEncodeCcreOverlaps requires non-empty seq_id and index_path"
                    .to_string(),
                cause_chain: vec![],
            });
        }
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
                cause_chain: vec![],
            })?;
        let sequence_length = dna.len();
        let start = start_0based.unwrap_or(0).min(sequence_length);
        let end = end_0based_exclusive
            .unwrap_or(sequence_length)
            .min(sequence_length);
        if end <= start {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "QueryEncodeCcreOverlaps requires a non-empty range (got {start}..{end})"
                ),
                cause_chain: vec![],
            });
        }
        let max_rows = limit.unwrap_or(DEFAULT_ENCODE_CCRE_QUERY_LIMIT);
        if max_rows == 0 || max_rows > MAX_ENCODE_CCRE_QUERY_LIMIT {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "ENCODE cCRE query limit must be within 1..={MAX_ENCODE_CCRE_QUERY_LIMIT}"
                ),
                cause_chain: vec![],
            });
        }
        let anchor = self.sequence_genome_anchor_summary(seq_id)?;
        let index =
            crate::encode_ccre::read_interval_index(index_path).map_err(|message| EngineError {
                code: ErrorCode::InvalidInput,
                message,
                cause_chain: vec![],
            })?;
        if !crate::encode_ccre::genome_id_matches_source(&anchor.genome_id, &index.source) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "ENCODE cCRE source '{}' is {} / {} (taxon {}), but sequence '{}' is anchored to genome '{}'; human and mouse resources must not be mixed",
                    index.source.source_id,
                    index.source.species_scientific_name,
                    index.source.assembly_name,
                    index.source.taxon_id,
                    seq_id,
                    anchor.genome_id
                ),
                cause_chain: vec![],
            });
        }
        let classes = Self::normalized_encode_ccre_classes(classes, &index.source.primary_classes)?;
        let resolved_bed_path =
            crate::encode_ccre::resolve_bed_path(index_path, &index, bed_path_override).map_err(
                |message| EngineError {
                    code: ErrorCode::InvalidInput,
                    message,
                    cause_chain: vec![],
                },
            )?;
        crate::encode_ccre::validate_bed_identity(&resolved_bed_path, &index).map_err(
            |message| EngineError {
                code: ErrorCode::InvalidInput,
                message,
                cause_chain: vec![],
            },
        )?;
        let (query_genomic_start, query_genomic_end) =
            Self::encode_ccre_genomic_query_span(&anchor, start, end)?;
        let intervals = crate::encode_ccre::query_overlaps(
            &index,
            &resolved_bed_path,
            &anchor.chromosome,
            query_genomic_start,
            query_genomic_end,
            &classes,
        )
        .map_err(|message| EngineError {
            code: ErrorCode::InvalidInput,
            message,
            cause_chain: vec![],
        })?;
        let mut rows = intervals
            .into_iter()
            .filter_map(|interval| {
                Self::project_encode_ccre_interval(
                    interval,
                    &anchor,
                    query_genomic_start,
                    query_genomic_end,
                )
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            left.local_start_0based
                .cmp(&right.local_start_0based)
                .then(
                    left.local_end_0based_exclusive
                        .cmp(&right.local_end_0based_exclusive),
                )
                .then(
                    left.interval
                        .ccre_accession
                        .cmp(&right.interval.ccre_accession),
                )
        });
        let matched_ccre_count = rows.len();
        let mut summaries = BTreeMap::<String, (usize, usize)>::new();
        for row in &rows {
            let summary = summaries
                .entry(row.interval.ccre_class.clone())
                .or_insert((0, 0));
            summary.0 = summary.0.saturating_add(1);
            summary.1 = summary.1.saturating_add(row.overlap_bp);
        }
        rows.truncate(max_rows);
        let index_sha256 = format!(
            "sha256:{}",
            crate::digest_utils::sha256_file_hex(Path::new(index_path)).map_err(|error| {
                EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not hash ENCODE cCRE index '{index_path}': {error}"),
                    cause_chain: vec![],
                }
            })?
        );
        let identity = format!(
            "{}\n{}\n{}\n{}\n{}:{}..{}:{}\n{}..{}\n{}\n{}",
            seq_id,
            index_sha256,
            index.bed_sha256,
            anchor.genome_id,
            anchor.chromosome,
            anchor.start_1based,
            anchor.end_1based,
            anchor.strand.unwrap_or('+'),
            start,
            end,
            classes.join(","),
            max_rows
        );
        Ok(EncodeCcreOverlapReport {
            schema: ENCODE_CCRE_OVERLAP_SCHEMA.to_string(),
            report_id: crate::digest_utils::short_sha256_id("encode_ccre_overlap", &identity),
            op_id: None,
            run_id: None,
            seq_id: seq_id.to_string(),
            index_path: index_path.to_string(),
            resolved_bed_path: resolved_bed_path.to_string_lossy().into_owned(),
            index_sha256,
            source: index.source,
            bed_sha256: index.bed_sha256,
            content_identity_verified: true,
            assembly_match_status: "assembly_and_species_matched".to_string(),
            genome_anchor: Some(EncodeCcreGenomeAnchor {
                seq_id: anchor.seq_id,
                genome_id: anchor.genome_id,
                chromosome: anchor.chromosome,
                start_1based: anchor.start_1based,
                end_1based: anchor.end_1based,
                strand: anchor.strand,
                anchor_verified: anchor.anchor_verified,
            }),
            query_start_0based: start,
            query_end_0based_exclusive: end,
            query_length_bp: end.saturating_sub(start),
            requested_classes: classes,
            max_rows,
            matched_ccre_count,
            returned_ccre_count: rows.len(),
            truncated: matched_ccre_count > rows.len(),
            class_summaries: summaries
                .into_iter()
                .map(|(ccre_class, (matched_count, overlap_bp))| EncodeCcreClassSummary {
                    ccre_class,
                    matched_count,
                    overlap_bp,
                })
                .collect(),
            rows,
            evidence_statement: ENCODE_CCRE_EVIDENCE_STATEMENT.to_string(),
            non_claims: vec![
                "A cCRE overlap does not prove enhancer activity in the experimental cell type."
                    .to_string(),
                "A cCRE overlap does not identify a regulated target gene or establish causal regulation."
                    .to_string(),
                "The ELS source excludes promoter-like and other Registry V4 classes."
                    .to_string(),
            ],
            warnings: vec![],
        })
    }

    fn is_generated_encode_ccre_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated")
            .any(|value| value.eq_ignore_ascii_case("encode_screen_ccre"))
    }

    fn encode_ccre_annotation_id(row: &EncodeCcreOverlapRow) -> String {
        format!(
            "{}:{}-{}:{}",
            row.interval.chromosome,
            row.interval.start_0based,
            row.interval.end_0based_exclusive,
            row.interval.ccre_accession
        )
    }

    fn build_encode_ccre_feature(
        row: &EncodeCcreOverlapRow,
        report: &EncodeCcreOverlapReport,
    ) -> gb_io::seq::Feature {
        let annotation_id = Self::encode_ccre_annotation_id(row);
        gb_io::seq::Feature {
            kind: "regulatory_region".into(),
            location: gb_io::seq::Location::simple_range(
                row.local_start_0based as i64,
                row.local_end_0based_exclusive as i64,
            ),
            qualifiers: vec![
                (
                    "label".into(),
                    Some(format!(
                        "{} ({})",
                        row.interval.ccre_accession, row.interval.ccre_class
                    )),
                ),
                (
                    "note".into(),
                    Some(
                        "ENCODE SCREEN Registry V4 candidate enhancer-like signature; overlap is evidence, not proof of enhancer activity or target-gene regulation."
                            .to_string(),
                    ),
                ),
                ("gentle_generated".into(), Some("encode_screen_ccre".to_string())),
                ("gentle_feature_source".into(), Some("ENCODE_SCREEN".to_string())),
                ("encode_ccre_annotation_id".into(), Some(annotation_id)),
                (
                    "encode_ccre_accession".into(),
                    Some(row.interval.ccre_accession.clone()),
                ),
                (
                    "encode_dhs_accession".into(),
                    Some(row.interval.dhs_accession.clone()),
                ),
                (
                    "encode_ccre_class".into(),
                    Some(row.interval.ccre_class.clone()),
                ),
                (
                    "regulatory_class".into(),
                    Some("enhancer_like_signature".to_string()),
                ),
                (
                    "encode_registry_version".into(),
                    Some(report.source.registry_version.clone()),
                ),
                (
                    "encode_source_id".into(),
                    Some(report.source.source_id.clone()),
                ),
                (
                    "encode_source_bed_sha256".into(),
                    Some(report.bed_sha256.clone()),
                ),
                (
                    "encode_genomic_chromosome".into(),
                    Some(row.interval.chromosome.clone()),
                ),
                (
                    "encode_genomic_start_0based".into(),
                    Some(row.interval.start_0based.to_string()),
                ),
                (
                    "encode_genomic_end_0based_exclusive".into(),
                    Some(row.interval.end_0based_exclusive.to_string()),
                ),
            ],
        }
    }

    pub fn materialize_encode_ccre_features(
        &mut self,
        seq_id: &str,
        index_path: &str,
        bed_path_override: Option<&str>,
        classes: &[String],
        max_features: Option<usize>,
        clear_existing: bool,
    ) -> Result<EncodeCcreMaterializationReport, EngineError> {
        let query = self.query_encode_ccre_overlaps(
            seq_id,
            index_path,
            bed_path_override,
            None,
            None,
            classes,
            max_features,
        )?;
        let matched_ccre_count = query.matched_ccre_count;
        let mut warnings = query.warnings.clone();
        if query.truncated {
            warnings.push(format!(
                "Materialized only the first {} of {} matching ENCODE cCRE intervals",
                query.returned_ccre_count, query.matched_ccre_count
            ));
        }
        let dna = self
            .state
            .sequences
            .get_mut(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
                cause_chain: vec![],
            })?;
        let mut removed_existing_count = 0;
        if clear_existing {
            let before = dna.features().len();
            dna.features_mut()
                .retain(|feature| !Self::is_generated_encode_ccre_feature(feature));
            removed_existing_count = before.saturating_sub(dna.features().len());
        }
        let mut existing_ids = dna
            .features()
            .iter()
            .filter(|feature| Self::is_generated_encode_ccre_feature(feature))
            .flat_map(|feature| feature.qualifier_values("encode_ccre_annotation_id"))
            .map(str::to_string)
            .collect::<HashSet<_>>();
        let mut skipped_existing_count = 0usize;
        let mut feature_ids = Vec::new();
        for row in &query.rows {
            let annotation_id = Self::encode_ccre_annotation_id(row);
            if !existing_ids.insert(annotation_id.clone()) {
                skipped_existing_count = skipped_existing_count.saturating_add(1);
                continue;
            }
            let feature_id = dna.features().len();
            dna.features_mut()
                .push(Self::build_encode_ccre_feature(row, &query));
            feature_ids.push(feature_id);
        }
        let added_feature_count = feature_ids.len();
        Self::prepare_sequence(dna);
        let report_id = crate::digest_utils::short_sha256_id(
            "encode_ccre_materialization",
            &format!(
                "{}\n{}\n{}\n{}\n{}\n{:?}",
                query.report_id,
                clear_existing,
                added_feature_count,
                skipped_existing_count,
                removed_existing_count,
                feature_ids,
            ),
        );
        Ok(EncodeCcreMaterializationReport {
            schema: ENCODE_CCRE_MATERIALIZATION_SCHEMA.to_string(),
            report_id,
            op_id: None,
            run_id: None,
            seq_id: seq_id.to_string(),
            index_path: index_path.to_string(),
            resolved_bed_path: query.resolved_bed_path,
            source: query.source,
            bed_sha256: query.bed_sha256,
            requested_classes: query.requested_classes,
            matched_ccre_count,
            added_feature_count,
            skipped_existing_count,
            removed_existing_count,
            feature_ids,
            evidence_statement: ENCODE_CCRE_EVIDENCE_STATEMENT.to_string(),
            warnings,
        })
    }
}
