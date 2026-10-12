//! Ensembl Regulation overlap projection and feature materialization.
//!
//! Assembly validation and evidence semantics stay in the engine so all
//! adapters consume the same annotation-only contract.

use super::*;

const DEFAULT_ENSEMBL_REGULATION_QUERY_LIMIT: usize = 10_000;
const MAX_ENSEMBL_REGULATION_QUERY_LIMIT: usize = 1_000_000;
const ENSEMBL_REGULATION_EVIDENCE_STATEMENT: &str = "Overlap with a release-bound Ensembl regulatory feature is regulatory-annotation evidence; activity and quantitative primary signal require separate evidence.";

impl GentleEngine {
    fn normalized_ensembl_regulation_types(
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
                    message: "Ensembl Regulation feature-type filters must not be blank"
                        .to_string(),
                    cause_chain: vec![],
                });
            }
            let Some(canonical) = available_by_lower.get(&token.to_ascii_lowercase()) else {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Ensembl Regulation feature type '{}' is not present in source; available types: {}",
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

    pub(crate) fn ensembl_regulation_genomic_query_span(
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
                message: "Ensembl Regulation query resolves to an empty genomic interval"
                    .to_string(),
                cause_chain: vec![],
            });
        }
        Ok(span)
    }

    pub(crate) fn project_ensembl_regulation_interval(
        interval: EnsemblRegulationInterval,
        anchor: &SequenceGenomeAnchorSummary,
        query_genomic_start: u64,
        query_genomic_end: u64,
    ) -> Option<EnsemblRegulationOverlapRow> {
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
        Some(EnsemblRegulationOverlapRow {
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

    #[allow(clippy::too_many_arguments)]
    pub fn query_ensembl_regulation_overlaps(
        &self,
        seq_id: &str,
        index_path: &str,
        intervals_path_override: Option<&str>,
        start_0based: Option<usize>,
        end_0based_exclusive: Option<usize>,
        feature_types: &[String],
        limit: Option<usize>,
    ) -> Result<EnsemblRegulationOverlapReport, EngineError> {
        let seq_id = seq_id.trim();
        if seq_id.is_empty() || index_path.trim().is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "QueryEnsemblRegulationOverlaps requires non-empty seq_id and index_path"
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
                    "QueryEnsemblRegulationOverlaps requires a non-empty range (got {start}..{end})"
                ),
                cause_chain: vec![],
            });
        }
        let max_rows = limit.unwrap_or(DEFAULT_ENSEMBL_REGULATION_QUERY_LIMIT);
        if max_rows == 0 || max_rows > MAX_ENSEMBL_REGULATION_QUERY_LIMIT {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Ensembl Regulation query limit must be within 1..={MAX_ENSEMBL_REGULATION_QUERY_LIMIT}"
                ),
                cause_chain: vec![],
            });
        }
        let anchor = self.sequence_genome_anchor_summary(seq_id)?;
        let index =
            crate::ensembl_regulation::read_interval_index(index_path).map_err(|message| {
                EngineError {
                    code: ErrorCode::InvalidInput,
                    message,
                    cause_chain: vec![],
                }
            })?;
        if !crate::ensembl_regulation::genome_id_matches_source(&anchor.genome_id, &index.source) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Ensembl Regulation source '{}' is {} / {} ({}, taxon {}), but sequence '{}' is anchored to genome '{}'; species and assemblies must not be mixed",
                    index.source.source_id,
                    index.source.species_scientific_name,
                    index.source.assembly_name,
                    index.source.assembly_accession,
                    index.source.taxon_id,
                    seq_id,
                    anchor.genome_id
                ),
                cause_chain: vec![],
            });
        }
        let feature_types =
            Self::normalized_ensembl_regulation_types(feature_types, &index.source.feature_types)?;
        let resolved_intervals_path = crate::ensembl_regulation::resolve_intervals_path(
            index_path,
            &index,
            intervals_path_override,
        )
        .map_err(|message| EngineError {
            code: ErrorCode::InvalidInput,
            message,
            cause_chain: vec![],
        })?;
        crate::ensembl_regulation::validate_intervals_identity(&resolved_intervals_path, &index)
            .map_err(|message| EngineError {
                code: ErrorCode::InvalidInput,
                message,
                cause_chain: vec![],
            })?;
        let (query_genomic_start, query_genomic_end) =
            Self::ensembl_regulation_genomic_query_span(&anchor, start, end)?;
        let intervals = crate::ensembl_regulation::query_overlaps(
            &index,
            &resolved_intervals_path,
            &anchor.chromosome,
            query_genomic_start,
            query_genomic_end,
            &feature_types,
        )
        .map_err(|message| EngineError {
            code: ErrorCode::InvalidInput,
            message,
            cause_chain: vec![],
        })?;
        let mut rows = intervals
            .into_iter()
            .filter_map(|interval| {
                Self::project_ensembl_regulation_interval(
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
                .then(left.interval.feature_id.cmp(&right.interval.feature_id))
        });
        let matched_feature_count = rows.len();
        let mut summaries = BTreeMap::<String, (usize, usize)>::new();
        for row in &rows {
            let summary = summaries
                .entry(row.interval.feature_type.clone())
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
                    message: format!(
                        "Could not hash Ensembl Regulation index '{index_path}': {error}"
                    ),
                    cause_chain: vec![],
                }
            })?
        );
        let identity = format!(
            "{}\n{}\n{}\n{}\n{}:{}..{}:{}\n{}..{}\n{}\n{}",
            seq_id,
            index_sha256,
            index.intervals_sha256,
            anchor.genome_id,
            anchor.chromosome,
            anchor.start_1based,
            anchor.end_1based,
            anchor.strand.unwrap_or('+'),
            start,
            end,
            feature_types.join(","),
            max_rows
        );
        Ok(EnsemblRegulationOverlapReport {
            schema: ENSEMBL_REGULATION_OVERLAP_SCHEMA.to_string(),
            report_id: crate::digest_utils::short_sha256_id(
                "ensembl_regulation_overlap",
                &identity,
            ),
            op_id: None,
            run_id: None,
            seq_id: seq_id.to_string(),
            index_path: index_path.to_string(),
            resolved_intervals_path: resolved_intervals_path.to_string_lossy().into_owned(),
            index_sha256,
            source: index.source,
            intervals_sha256: index.intervals_sha256,
            content_identity_verified: true,
            assembly_match_status: "assembly_and_species_matched_source_release_pinned"
                .to_string(),
            genome_anchor: Some(EnsemblRegulationGenomeAnchor {
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
            requested_feature_types: feature_types,
            max_rows,
            matched_feature_count,
            returned_feature_count: rows.len(),
            truncated: matched_feature_count > rows.len(),
            type_summaries: summaries
                .into_iter()
                .map(
                    |(feature_type, (matched_count, overlap_bp))| EnsemblRegulationTypeSummary {
                        feature_type,
                        matched_count,
                        overlap_bp,
                    },
                )
                .collect(),
            rows,
            evidence_statement: ENSEMBL_REGULATION_EVIDENCE_STATEMENT.to_string(),
            non_claims: vec![
                "Regulatory-feature overlap does not establish activity in the experimental biosample."
                    .to_string(),
                "An associated-gene annotation does not by itself establish causal regulation."
                    .to_string(),
                "Regulatory activity and epigenome-specific quantitative signals were not evaluated by this annotation query."
                    .to_string(),
            ],
            warnings: vec![],
        })
    }

    fn is_generated_ensembl_regulation_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated")
            .any(|value| value.eq_ignore_ascii_case("ensembl_regulation"))
    }

    fn ensembl_regulation_annotation_id(row: &EnsemblRegulationOverlapRow) -> String {
        format!(
            "{}:{}-{}:{}",
            row.interval.chromosome,
            row.interval.start_0based,
            row.interval.end_0based_exclusive,
            row.interval.feature_id
        )
    }

    fn build_ensembl_regulation_feature(
        row: &EnsemblRegulationOverlapRow,
        report: &EnsemblRegulationOverlapReport,
    ) -> gb_io::seq::Feature {
        let annotation_id = Self::ensembl_regulation_annotation_id(row);
        let mut qualifiers = vec![
            (
                "label".into(),
                Some(format!(
                    "{} ({})",
                    row.interval.feature_id, row.interval.feature_type
                )),
            ),
            (
                "note".into(),
                Some(
                    "Ensembl regulatory annotation; overlap is evidence, while activity and quantitative primary signal require separate evidence."
                        .to_string(),
                ),
            ),
            (
                "gentle_generated".into(),
                Some("ensembl_regulation".to_string()),
            ),
            (
                "gentle_feature_source".into(),
                Some("Ensembl_Regulation".to_string()),
            ),
            (
                "ensembl_regulation_annotation_id".into(),
                Some(annotation_id),
            ),
            (
                "ensembl_regulatory_feature_id".into(),
                Some(row.interval.feature_id.clone()),
            ),
            (
                "regulatory_class".into(),
                Some(row.interval.feature_type.clone()),
            ),
            (
                "ensembl_regulation_release".into(),
                Some(report.source.annotation_release.clone()),
            ),
            (
                "ensembl_regulation_source_id".into(),
                Some(report.source.source_id.clone()),
            ),
            (
                "ensembl_regulation_intervals_sha256".into(),
                Some(report.intervals_sha256.clone()),
            ),
            (
                "ensembl_genomic_chromosome".into(),
                Some(row.interval.chromosome.clone()),
            ),
            (
                "ensembl_genomic_start_0based".into(),
                Some(row.interval.start_0based.to_string()),
            ),
            (
                "ensembl_genomic_end_0based_exclusive".into(),
                Some(row.interval.end_0based_exclusive.to_string()),
            ),
        ];
        if let Some(value) = row.interval.extended_start_0based {
            qualifiers.push((
                "ensembl_extended_start_0based".into(),
                Some(value.to_string()),
            ));
        }
        if let Some(value) = row.interval.extended_end_0based_exclusive {
            qualifiers.push((
                "ensembl_extended_end_0based_exclusive".into(),
                Some(value.to_string()),
            ));
        }
        qualifiers.extend(
            row.interval
                .associated_gene_ids
                .iter()
                .cloned()
                .map(|value| ("ensembl_associated_gene_id".into(), Some(value))),
        );
        qualifiers.extend(
            row.interval
                .associated_gene_names
                .iter()
                .cloned()
                .map(|value| ("ensembl_associated_gene_name".into(), Some(value))),
        );
        if let Some(strand) = row.interval.strand {
            qualifiers.push(("ensembl_genomic_strand".into(), Some(strand.to_string())));
        }
        let location = gb_io::seq::Location::simple_range(
            row.local_start_0based as i64,
            row.local_end_0based_exclusive as i64,
        );
        let anchor_is_reverse = report
            .genome_anchor
            .as_ref()
            .is_some_and(|anchor| anchor.strand == Some('-'));
        let local_is_reverse = row
            .interval
            .strand
            .is_some_and(|strand| (strand == '-') != anchor_is_reverse);
        gb_io::seq::Feature {
            kind: "regulatory_region".into(),
            location: if local_is_reverse {
                gb_io::seq::Location::Complement(Box::new(location))
            } else {
                location
            },
            qualifiers,
        }
    }

    pub fn materialize_ensembl_regulation_features(
        &mut self,
        seq_id: &str,
        index_path: &str,
        intervals_path_override: Option<&str>,
        feature_types: &[String],
        max_features: Option<usize>,
        clear_existing: bool,
    ) -> Result<EnsemblRegulationMaterializationReport, EngineError> {
        let query = self.query_ensembl_regulation_overlaps(
            seq_id,
            index_path,
            intervals_path_override,
            None,
            None,
            feature_types,
            max_features,
        )?;
        let matched_feature_count = query.matched_feature_count;
        let mut warnings = query.warnings.clone();
        if query.truncated {
            warnings.push(format!(
                "Materialized only the first {} of {} matching Ensembl regulatory features",
                query.returned_feature_count, query.matched_feature_count
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
                .retain(|feature| !Self::is_generated_ensembl_regulation_feature(feature));
            removed_existing_count = before.saturating_sub(dna.features().len());
        }
        let mut existing_ids = dna
            .features()
            .iter()
            .filter(|feature| Self::is_generated_ensembl_regulation_feature(feature))
            .flat_map(|feature| feature.qualifier_values("ensembl_regulation_annotation_id"))
            .map(str::to_string)
            .collect::<HashSet<_>>();
        let mut skipped_existing_count = 0usize;
        let mut feature_ids = Vec::new();
        for row in &query.rows {
            let annotation_id = Self::ensembl_regulation_annotation_id(row);
            if !existing_ids.insert(annotation_id) {
                skipped_existing_count = skipped_existing_count.saturating_add(1);
                continue;
            }
            let feature_id = dna.features().len();
            dna.features_mut()
                .push(Self::build_ensembl_regulation_feature(row, &query));
            feature_ids.push(feature_id);
        }
        let added_feature_count = feature_ids.len();
        Self::prepare_sequence(dna);
        let report_id = crate::digest_utils::short_sha256_id(
            "ensembl_regulation_materialization",
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
        Ok(EnsemblRegulationMaterializationReport {
            schema: ENSEMBL_REGULATION_MATERIALIZATION_SCHEMA.to_string(),
            report_id,
            op_id: None,
            run_id: None,
            seq_id: seq_id.to_string(),
            index_path: index_path.to_string(),
            resolved_intervals_path: query.resolved_intervals_path,
            source: query.source,
            intervals_sha256: query.intervals_sha256,
            requested_feature_types: query.requested_feature_types,
            matched_feature_count,
            added_feature_count,
            skipped_existing_count,
            removed_existing_count,
            feature_ids,
            evidence_statement: ENSEMBL_REGULATION_EVIDENCE_STATEMENT.to_string(),
            warnings,
        })
    }
}
