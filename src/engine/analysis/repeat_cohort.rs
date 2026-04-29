//! RepeatMasker-backed repeat cohorts and transcript-aware window geometry.
//!
//! This module keeps UCSC `rmsk` parsing and repeat-family cohort construction
//! in the shared engine so GUI, CLI, and agent-facing adapters can reuse one
//! deterministic interpretation path.

use super::*;

impl GentleEngine {
    fn open_rmsk_reader(path: &str) -> Result<Box<dyn BufRead>, EngineError> {
        let file = File::open(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not open RepeatMasker file '{}': {e}", path),
        })?;
        if path.to_ascii_lowercase().ends_with(".gz") {
            Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
        } else {
            Ok(Box::new(BufReader::new(file)))
        }
    }

    fn normalize_repeat_strand(raw: &str) -> String {
        match raw.trim() {
            "-" | "C" | "c" => "-".to_string(),
            "+" => "+".to_string(),
            other if other.is_empty() => ".".to_string(),
            other => other.to_string(),
        }
    }

    pub(crate) fn normalized_repeat_alias(
        rep_name: &str,
        rep_class: &str,
        rep_family: &str,
    ) -> String {
        crate::ucsc_rmsk::normalized_repeat_alias(rep_name, rep_class, rep_family)
    }

    fn parse_optional_f64(raw: Option<&str>) -> Option<f64> {
        raw.and_then(|value| value.trim().parse::<f64>().ok())
    }

    fn parse_ucsc_rmsk_record(
        genome_id: &str,
        line: &str,
        line_number: usize,
    ) -> Result<Option<RepeatAnnotationRecord>, String> {
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
        {
            return Ok(None);
        }
        let fields = trimmed.split_whitespace().collect::<Vec<_>>();
        if fields.is_empty() {
            return Ok(None);
        }
        let first = fields[0].to_ascii_lowercase();
        if matches!(first.as_str(), "sw" | "swscore" | "bin") {
            return Ok(None);
        }

        let parsed = if fields.len() >= 17 {
            let chromosome = fields[5].to_string();
            let start_0based = fields[6]
                .parse::<usize>()
                .map_err(|_| format!("invalid genoStart '{}'", fields[6]))?;
            let end_0based_exclusive = fields[7]
                .parse::<usize>()
                .map_err(|_| format!("invalid genoEnd '{}'", fields[7]))?;
            (
                chromosome,
                start_0based,
                end_0based_exclusive,
                Self::normalize_repeat_strand(fields[9]),
                fields[10].to_string(),
                fields[11].to_string(),
                fields[12].to_string(),
                Self::parse_optional_f64(fields.get(1).copied()),
                Self::parse_optional_f64(fields.get(2).copied()),
            )
        } else if fields.len() >= 16 && fields[0].parse::<f64>().is_ok() {
            let chromosome = fields[4].to_string();
            let start_0based = fields[5]
                .parse::<usize>()
                .map_err(|_| format!("invalid genoStart '{}'", fields[5]))?;
            let end_0based_exclusive = fields[6]
                .parse::<usize>()
                .map_err(|_| format!("invalid genoEnd '{}'", fields[6]))?;
            (
                chromosome,
                start_0based,
                end_0based_exclusive,
                Self::normalize_repeat_strand(fields[8]),
                fields[9].to_string(),
                fields[10].to_string(),
                fields[11].to_string(),
                Self::parse_optional_f64(fields.first().copied()),
                Self::parse_optional_f64(fields.get(1).copied()),
            )
        } else if fields.len() >= 8 {
            let chromosome = fields[0].to_string();
            let start_0based = fields[1]
                .parse::<usize>()
                .map_err(|_| format!("invalid chromStart '{}'", fields[1]))?;
            let end_0based_exclusive = fields[2]
                .parse::<usize>()
                .map_err(|_| format!("invalid chromEnd '{}'", fields[2]))?;
            (
                chromosome,
                start_0based,
                end_0based_exclusive,
                Self::normalize_repeat_strand(fields[5]),
                fields[3].to_string(),
                fields[6].to_string(),
                fields[7].to_string(),
                Self::parse_optional_f64(fields.get(4).copied()),
                None,
            )
        } else {
            return Err("expected UCSC rmsk table or BED6+repClass+repFamily row".to_string());
        };

        let (
            chromosome,
            start_0based,
            end_0based_exclusive,
            strand,
            rep_name,
            rep_class,
            rep_family,
            score,
            milli_div,
        ) = parsed;
        if end_0based_exclusive <= start_0based {
            return Err(format!(
                "end {} must be greater than start {}",
                end_0based_exclusive, start_0based
            ));
        }
        let normalized_alias = Self::normalized_repeat_alias(&rep_name, &rep_class, &rep_family);
        Ok(Some(RepeatAnnotationRecord {
            annotation_id: format!(
                "{}:{}-{}:{}:{}",
                chromosome, start_0based, end_0based_exclusive, rep_name, line_number
            ),
            genome_id: genome_id.to_string(),
            chromosome,
            start_0based,
            end_0based_exclusive,
            start_1based: start_0based.saturating_add(1),
            end_1based: end_0based_exclusive,
            strand,
            rep_name,
            rep_class,
            rep_family,
            normalized_alias,
            score,
            milli_div,
            source_line_number: Some(line_number),
        }))
    }

    fn repeat_filter_token_matches(tokens: &[String], value: &str) -> bool {
        tokens.is_empty()
            || tokens
                .iter()
                .map(|token| token.trim())
                .filter(|token| !token.is_empty())
                .any(|token| token.eq_ignore_ascii_case(value))
    }

    fn repeat_annotation_matches_filter(
        row: &RepeatAnnotationRecord,
        filter: &RepeatAnnotationFilter,
    ) -> bool {
        if !Self::repeat_filter_token_matches(&filter.rep_names, &row.rep_name)
            || !Self::repeat_filter_token_matches(&filter.rep_classes, &row.rep_class)
            || !Self::repeat_filter_token_matches(&filter.rep_families, &row.rep_family)
            || !Self::repeat_filter_token_matches(&filter.normalized_aliases, &row.normalized_alias)
        {
            return false;
        }
        if let Some(chromosome) = filter.chromosome.as_deref() {
            if !chromosome.trim().is_empty() && !chromosome.eq_ignore_ascii_case(&row.chromosome) {
                return false;
            }
        }
        if filter.span_start_0based.is_some() || filter.span_end_0based_exclusive.is_some() {
            let span_start = filter.span_start_0based.unwrap_or(0);
            let span_end = filter.span_end_0based_exclusive.unwrap_or(usize::MAX);
            if row.start_0based >= span_end || row.end_0based_exclusive <= span_start {
                return false;
            }
        }
        true
    }

    pub(crate) fn query_repeat_annotations(
        &self,
        genome_id: &str,
        rmsk_path: &str,
        filter: &RepeatAnnotationFilter,
        limit: Option<usize>,
    ) -> Result<RepeatAnnotationQueryReport, EngineError> {
        if genome_id.trim().is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "QueryRepeatAnnotations requires non-empty genome_id".to_string(),
            });
        }
        if rmsk_path.trim().is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "QueryRepeatAnnotations requires non-empty rmsk_path".to_string(),
            });
        }
        let reader = Self::open_rmsk_reader(rmsk_path)?;
        let mut parsed_row_count = 0usize;
        let mut matched_row_count = 0usize;
        let mut malformed_line_count = 0usize;
        let mut malformed_examples = vec![];
        let mut rows = vec![];
        let effective_limit = limit.unwrap_or(usize::MAX);
        for (idx, line) in reader.lines().enumerate() {
            let line_number = idx + 1;
            let line = line.map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not read RepeatMasker line {line_number}: {e}"),
            })?;
            match Self::parse_ucsc_rmsk_record(genome_id, &line, line_number) {
                Ok(Some(row)) => {
                    parsed_row_count = parsed_row_count.saturating_add(1);
                    if Self::repeat_annotation_matches_filter(&row, filter) {
                        matched_row_count = matched_row_count.saturating_add(1);
                        if rows.len() < effective_limit {
                            rows.push(row);
                        }
                    }
                }
                Ok(None) => {}
                Err(err) => {
                    malformed_line_count = malformed_line_count.saturating_add(1);
                    if malformed_examples.len() < 5 {
                        malformed_examples.push(format!("line {line_number}: {err}"));
                    }
                }
            }
        }
        Ok(RepeatAnnotationQueryReport {
            schema: REPEAT_ANNOTATION_QUERY_REPORT_SCHEMA.to_string(),
            genome_id: genome_id.to_string(),
            rmsk_path: rmsk_path.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            filter: filter.clone(),
            parsed_row_count,
            matched_row_count,
            returned_row_count: rows.len(),
            malformed_line_count,
            malformed_examples,
            rows,
            warnings: vec![],
        })
    }

    fn repeat_record_from_indexed_interval(
        genome_id: &str,
        interval: &crate::ucsc_rmsk::UcscRmskIndexedInterval,
    ) -> RepeatAnnotationRecord {
        RepeatAnnotationRecord {
            annotation_id: format!(
                "{}:{}-{}:{}:{}",
                interval.chromosome,
                interval.start_0based,
                interval.end_0based_exclusive,
                interval.rep_name,
                interval.source_row + 1
            ),
            genome_id: genome_id.to_string(),
            chromosome: interval.chromosome.clone(),
            start_0based: interval.start_0based,
            end_0based_exclusive: interval.end_0based_exclusive,
            start_1based: interval.start_0based.saturating_add(1),
            end_1based: interval.end_0based_exclusive,
            strand: interval.strand.clone(),
            rep_name: interval.rep_name.clone(),
            rep_class: interval.rep_class.clone(),
            rep_family: interval.rep_family.clone(),
            normalized_alias: interval.normalized_alias.clone(),
            score: Some(interval.sw_score as f64),
            milli_div: Some(interval.milli_div as f64),
            source_line_number: Some(interval.source_row + 1),
        }
    }

    fn rmsk_projection_anchor_from_sequence(
        anchor: &SequenceGenomeAnchorSummary,
    ) -> crate::ucsc_rmsk::UcscRmskProjectionAnchor {
        crate::ucsc_rmsk::UcscRmskProjectionAnchor {
            chromosome: anchor.chromosome.clone(),
            start_1based: anchor.start_1based,
            end_1based: anchor.end_1based,
            strand: anchor.strand,
        }
    }

    pub fn query_sequence_repeat_overlaps(
        &self,
        seq_id: &str,
        rmsk_index_path: &str,
        start_0based: Option<usize>,
        end_0based_exclusive: Option<usize>,
        limit: Option<usize>,
    ) -> Result<SequenceRepeatOverlapReport, EngineError> {
        let seq_id = seq_id.trim();
        if seq_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "QueryRepeatOverlaps requires non-empty seq_id".to_string(),
            });
        }
        if rmsk_index_path.trim().is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "QueryRepeatOverlaps requires non-empty rmsk_index_path".to_string(),
            });
        }
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let sequence_length = dna.len();
        let start = start_0based.unwrap_or(0).min(sequence_length);
        let end = end_0based_exclusive
            .unwrap_or(sequence_length)
            .min(sequence_length);
        if end < start {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "QueryRepeatOverlaps requires end >= start (got {}..{})",
                    start, end
                ),
            });
        }
        let anchor = self.sequence_genome_anchor_summary(seq_id)?;
        let index =
            crate::ucsc_rmsk::read_ucsc_rmsk_interval_index(rmsk_index_path).map_err(|e| {
                EngineError {
                    code: ErrorCode::InvalidInput,
                    message: e,
                }
            })?;
        let projection_anchor = Self::rmsk_projection_anchor_from_sequence(&anchor);
        let overlaps = crate::ucsc_rmsk::project_index_overlaps_to_anchor(
            &index,
            &projection_anchor,
            Some(start),
            Some(end),
            None,
        );
        let matched_repeat_count = overlaps.len();
        let effective_limit = limit.unwrap_or(usize::MAX);
        let rows = overlaps
            .into_iter()
            .take(effective_limit)
            .map(|overlap| SequenceRepeatOverlapRow {
                repeat: Self::repeat_record_from_indexed_interval(
                    &anchor.genome_id,
                    &overlap.interval,
                ),
                local_start_0based: overlap.local_start_0based,
                local_end_0based_exclusive: overlap.local_end_0based_exclusive,
                local_strand: overlap.local_strand,
                genomic_start_0based: overlap.genomic_start_0based,
                genomic_end_0based_exclusive: overlap.genomic_end_0based_exclusive,
                overlap_bp: overlap.overlap_bp,
                clipped: overlap.clipped,
            })
            .collect::<Vec<_>>();
        Ok(SequenceRepeatOverlapReport {
            schema: SEQUENCE_REPEAT_OVERLAP_REPORT_SCHEMA.to_string(),
            seq_id: seq_id.to_string(),
            rmsk_index_path: rmsk_index_path.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            genome_anchor: Some(anchor),
            query_start_0based: Some(start),
            query_end_0based_exclusive: Some(end),
            matched_repeat_count,
            returned_repeat_count: rows.len(),
            rows,
            warnings: vec![],
        })
    }

    fn is_generated_ucsc_rmsk_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated")
            .any(|value| value.eq_ignore_ascii_case("ucsc_rmsk"))
    }

    fn repeat_feature_annotation_id(
        interval: &crate::ucsc_rmsk::UcscRmskIndexedInterval,
    ) -> String {
        format!(
            "{}:{}-{}:{}:{}",
            interval.chromosome,
            interval.start_0based,
            interval.end_0based_exclusive,
            interval.rep_name,
            interval.source_row + 1
        )
    }

    fn build_ucsc_rmsk_feature(
        overlap: &crate::ucsc_rmsk::UcscRmskProjectedOverlap,
        rmsk_index_path: &str,
    ) -> gb_io::seq::Feature {
        let mut location = gb_io::seq::Location::simple_range(
            overlap.local_start_0based as i64,
            overlap.local_end_0based_exclusive as i64,
        );
        if overlap.local_strand == "-" {
            location = gb_io::seq::Location::Complement(Box::new(location));
        }
        let interval = &overlap.interval;
        let annotation_id = Self::repeat_feature_annotation_id(interval);
        let label = if interval.rep_family.trim().is_empty() {
            interval.rep_name.clone()
        } else {
            format!("{} ({})", interval.rep_name, interval.rep_family)
        };
        gb_io::seq::Feature {
            kind: "repeat_region".into(),
            location,
            qualifiers: vec![
                ("label".into(), Some(label)),
                (
                    "note".into(),
                    Some("UCSC RepeatMasker rmsk interval".to_string()),
                ),
                ("gentle_generated".into(), Some("ucsc_rmsk".to_string())),
                (
                    "gentle_feature_source".into(),
                    Some("ucsc_rmsk".to_string()),
                ),
                ("rmsk_index_path".into(), Some(rmsk_index_path.to_string())),
                ("rmsk_annotation_id".into(), Some(annotation_id)),
                (
                    "rmsk_row_offset".into(),
                    Some(interval.source_row.to_string()),
                ),
                ("rmsk_bin".into(), Some(interval.bin.to_string())),
                ("rmsk_name".into(), Some(interval.rep_name.clone())),
                ("rmsk_class".into(), Some(interval.rep_class.clone())),
                ("rmsk_family".into(), Some(interval.rep_family.clone())),
                ("repName".into(), Some(interval.rep_name.clone())),
                ("repClass".into(), Some(interval.rep_class.clone())),
                ("repFamily".into(), Some(interval.rep_family.clone())),
                ("repeat_name".into(), Some(interval.rep_name.clone())),
                ("repeat_class".into(), Some(interval.rep_class.clone())),
                ("repeat_family".into(), Some(interval.rep_family.clone())),
                (
                    "repeat_alias".into(),
                    Some(interval.normalized_alias.clone()),
                ),
                ("strand".into(), Some(overlap.local_strand.clone())),
                ("rmsk_genomic_strand".into(), Some(interval.strand.clone())),
                ("chromosome".into(), Some(interval.chromosome.clone())),
                (
                    "genomic_start_1based".into(),
                    Some(overlap.genomic_start_0based.saturating_add(1).to_string()),
                ),
                (
                    "genomic_end_1based".into(),
                    Some(overlap.genomic_end_0based_exclusive.to_string()),
                ),
                ("swScore".into(), Some(interval.sw_score.to_string())),
                ("milliDiv".into(), Some(interval.milli_div.to_string())),
                ("milliDel".into(), Some(interval.milli_del.to_string())),
                ("milliIns".into(), Some(interval.milli_ins.to_string())),
                ("rmsk_clipped".into(), Some(overlap.clipped.to_string())),
            ],
        }
    }

    pub fn materialize_sequence_repeat_features(
        &mut self,
        seq_id: &str,
        rmsk_index_path: &str,
        max_features: Option<usize>,
        clear_existing: bool,
    ) -> Result<RepeatFeatureMaterializationReport, EngineError> {
        let anchor = self.sequence_genome_anchor_summary(seq_id)?;
        let projection_anchor = Self::rmsk_projection_anchor_from_sequence(&anchor);
        let index =
            crate::ucsc_rmsk::read_ucsc_rmsk_interval_index(rmsk_index_path).map_err(|e| {
                EngineError {
                    code: ErrorCode::InvalidInput,
                    message: e,
                }
            })?;
        let overlaps = crate::ucsc_rmsk::project_index_overlaps_to_anchor(
            &index,
            &projection_anchor,
            None,
            None,
            None,
        );
        let mut warnings = vec![];
        let matched_repeat_count = overlaps.len();
        let feature_limit = max_features.unwrap_or(usize::MAX);
        if matched_repeat_count > feature_limit {
            warnings.push(format!(
                "Materialized only the first {} of {} overlapping UCSC rmsk repeat(s)",
                feature_limit, matched_repeat_count
            ));
        }
        let dna = self
            .state
            .sequences
            .get_mut(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let mut removed_existing_count = 0usize;
        if clear_existing {
            let before = dna.features().len();
            dna.features_mut()
                .retain(|feature| !Self::is_generated_ucsc_rmsk_feature(feature));
            removed_existing_count = before.saturating_sub(dna.features().len());
        }
        let existing_ids = dna
            .features()
            .iter()
            .filter(|feature| Self::is_generated_ucsc_rmsk_feature(feature))
            .flat_map(|feature| feature.qualifier_values("rmsk_annotation_id"))
            .map(str::to_string)
            .collect::<HashSet<_>>();
        let sequence_length = dna.len();
        let mut feature_ids = vec![];
        let mut skipped_existing_count = 0usize;
        for row in overlaps.into_iter().take(feature_limit) {
            let annotation_id = Self::repeat_feature_annotation_id(&row.interval);
            if row.local_start_0based >= sequence_length
                || row.local_end_0based_exclusive > sequence_length
                || row.local_end_0based_exclusive <= row.local_start_0based
            {
                warnings.push(format!(
                    "Skipped repeat '{}' because projected local span {}..{} is outside sequence length {}",
                    annotation_id,
                    row.local_start_0based,
                    row.local_end_0based_exclusive,
                    sequence_length
                ));
                continue;
            }
            if !clear_existing && existing_ids.contains(&annotation_id) {
                skipped_existing_count = skipped_existing_count.saturating_add(1);
                continue;
            }
            let feature_id = dna.features().len();
            dna.features_mut()
                .push(Self::build_ucsc_rmsk_feature(&row, rmsk_index_path));
            feature_ids.push(feature_id);
        }
        let added_feature_count = feature_ids.len();
        Self::prepare_sequence(dna);
        Ok(RepeatFeatureMaterializationReport {
            schema: REPEAT_FEATURE_MATERIALIZATION_REPORT_SCHEMA.to_string(),
            seq_id: seq_id.to_string(),
            rmsk_index_path: rmsk_index_path.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            genome_anchor: Some(anchor),
            matched_repeat_count,
            added_feature_count,
            skipped_existing_count,
            removed_existing_count,
            feature_ids,
            warnings,
        })
    }

    fn repeat_midpoint_1based(row: &RepeatAnnotationRecord) -> usize {
        row.start_1based
            .saturating_add(row.end_1based)
            .saturating_div(2)
    }

    fn oriented_window_around_anchor(
        chromosome: &str,
        anchor_1based: usize,
        strand: &str,
        upstream_bp: usize,
        downstream_bp: usize,
        mode: RepeatEnvironmentGeometryMode,
        reason: &str,
        transcript_id: Option<String>,
    ) -> RepeatEnvironmentGeometryWindow {
        let (start, end) = if strand == "-" {
            (
                anchor_1based.saturating_sub(downstream_bp).max(1),
                anchor_1based.saturating_add(upstream_bp),
            )
        } else {
            (
                anchor_1based.saturating_sub(upstream_bp).max(1),
                anchor_1based.saturating_add(downstream_bp),
            )
        };
        RepeatEnvironmentGeometryWindow {
            mode,
            available: true,
            reason: reason.to_string(),
            chromosome: Some(chromosome.to_string()),
            start_1based: Some(start),
            end_1based: Some(end.max(start)),
            anchor_1based: Some(anchor_1based),
            strand: Some(strand.to_string()),
            transcript_id,
        }
    }

    fn unavailable_repeat_geometry(
        mode: RepeatEnvironmentGeometryMode,
        reason: &str,
    ) -> RepeatEnvironmentGeometryWindow {
        RepeatEnvironmentGeometryWindow {
            mode,
            available: false,
            reason: reason.to_string(),
            ..RepeatEnvironmentGeometryWindow::default()
        }
    }

    fn signed_distance_on_strand(
        position_1based: usize,
        anchor_1based: usize,
        strand: &str,
    ) -> i64 {
        if strand == "-" {
            anchor_1based as i64 - position_1based as i64
        } else {
            position_1based as i64 - anchor_1based as i64
        }
    }

    pub(crate) fn transcript_context_for_repeat(
        repeat: &RepeatAnnotationRecord,
        transcript: &GenomeTranscriptRecord,
    ) -> RepeatTranscriptContext {
        let strand = transcript.strand.unwrap_or('+').to_string();
        let tss_1based = if strand == "-" {
            transcript.transcript_end_1based
        } else {
            transcript.transcript_start_1based
        };
        let repeat_midpoint = Self::repeat_midpoint_1based(repeat);
        let overlaps_repeat = transcript
            .chromosome
            .eq_ignore_ascii_case(&repeat.chromosome)
            && transcript.transcript_start_1based <= repeat.end_1based
            && transcript.transcript_end_1based >= repeat.start_1based;
        let repeat_relation = if overlaps_repeat {
            "overlaps_transcript".to_string()
        } else if Self::signed_distance_on_strand(repeat_midpoint, tss_1based, &strand) < 0 {
            "upstream_of_tss".to_string()
        } else {
            "downstream_of_tss".to_string()
        };
        let cds_5prime_boundary = if transcript.cds_1based.is_empty() {
            None
        } else if strand == "-" {
            transcript.cds_1based.iter().map(|(_, end)| *end).max()
        } else {
            transcript.cds_1based.iter().map(|(start, _)| *start).min()
        };
        let cds_stop_1based = if transcript.cds_1based.is_empty() {
            None
        } else if strand == "-" {
            transcript.cds_1based.iter().map(|(start, _)| *start).min()
        } else {
            transcript.cds_1based.iter().map(|(_, end)| *end).max()
        };
        let (inferred_5utr_start_1based, inferred_5utr_end_1based) =
            match (cds_5prime_boundary, strand.as_str()) {
                (Some(cds), "+") if transcript.transcript_start_1based < cds => (
                    Some(transcript.transcript_start_1based),
                    Some(cds.saturating_sub(1)),
                ),
                (Some(cds), "-") if cds < transcript.transcript_end_1based => (
                    Some(cds.saturating_add(1)),
                    Some(transcript.transcript_end_1based),
                ),
                _ => (None, None),
            };
        RepeatTranscriptContext {
            gene_label: transcript
                .gene_name
                .clone()
                .or_else(|| transcript.gene_id.clone())
                .unwrap_or_else(|| transcript.transcript_id.clone()),
            gene_id: transcript.gene_id.clone(),
            transcript_id: transcript.transcript_id.clone(),
            chromosome: transcript.chromosome.clone(),
            strand: strand.clone(),
            transcript_start_1based: transcript.transcript_start_1based,
            transcript_end_1based: transcript.transcript_end_1based,
            overlaps_repeat,
            repeat_relation,
            tss_1based: Some(tss_1based),
            inferred_5utr_start_1based,
            inferred_5utr_end_1based,
            cds_stop_1based,
            signed_tss_distance_bp: Some(Self::signed_distance_on_strand(
                repeat_midpoint,
                tss_1based,
                &strand,
            )),
            signed_stop_distance_bp: cds_stop_1based
                .map(|stop| Self::signed_distance_on_strand(repeat_midpoint, stop, &strand)),
        }
    }

    fn choose_primary_repeat_transcript_context(
        contexts: &[RepeatTranscriptContext],
    ) -> Option<&RepeatTranscriptContext> {
        contexts.iter().min_by(|left, right| {
            let left_overlap_rank = if left.overlaps_repeat { 0 } else { 1 };
            let right_overlap_rank = if right.overlaps_repeat { 0 } else { 1 };
            left_overlap_rank
                .cmp(&right_overlap_rank)
                .then_with(|| {
                    left.signed_tss_distance_bp
                        .unwrap_or(i64::MAX)
                        .abs()
                        .cmp(&right.signed_tss_distance_bp.unwrap_or(i64::MAX).abs())
                })
                .then(left.transcript_id.cmp(&right.transcript_id))
        })
    }

    pub(crate) fn build_repeat_geometry_windows(
        repeat: &RepeatAnnotationRecord,
        contexts: &[RepeatTranscriptContext],
        upstream_bp: usize,
        downstream_bp: usize,
    ) -> Vec<RepeatEnvironmentGeometryWindow> {
        let repeat_strand = if repeat.strand == "-" { "-" } else { "+" };
        let repeat_anchor = Self::repeat_midpoint_1based(repeat);
        let mut windows = vec![Self::oriented_window_around_anchor(
            &repeat.chromosome,
            repeat_anchor,
            repeat_strand,
            upstream_bp,
            downstream_bp,
            RepeatEnvironmentGeometryMode::RepeatMidpoint,
            "repeat midpoint from rmsk span",
            None,
        )];
        let primary = Self::choose_primary_repeat_transcript_context(contexts);
        if let Some(context) = primary {
            let tss = context.tss_1based;
            windows.push(if let Some(anchor) = tss {
                let reason = if context.inferred_5utr_start_1based.is_some() {
                    "transcript 5' boundary with inferred 5'UTR from CDS"
                } else {
                    "transcript 5' boundary; CDS/UTR inference unavailable"
                };
                Self::oriented_window_around_anchor(
                    &context.chromosome,
                    anchor,
                    &context.strand,
                    upstream_bp,
                    downstream_bp,
                    RepeatEnvironmentGeometryMode::Transcript5utrStart,
                    reason,
                    Some(context.transcript_id.clone()),
                )
            } else {
                Self::unavailable_repeat_geometry(
                    RepeatEnvironmentGeometryMode::Transcript5utrStart,
                    "no transcript TSS context available",
                )
            });
            windows.push(if let Some(anchor) = tss {
                Self::oriented_window_around_anchor(
                    &context.chromosome,
                    anchor,
                    &context.strand,
                    upstream_bp,
                    downstream_bp,
                    RepeatEnvironmentGeometryMode::Pol2PromoterUpstream,
                    "strand-aware promoter-side window upstream of transcript 5' boundary",
                    Some(context.transcript_id.clone()),
                )
            } else {
                Self::unavailable_repeat_geometry(
                    RepeatEnvironmentGeometryMode::Pol2PromoterUpstream,
                    "no transcript TSS context available",
                )
            });
            windows.push(if let Some(anchor) = context.cds_stop_1based {
                Self::oriented_window_around_anchor(
                    &context.chromosome,
                    anchor,
                    &context.strand,
                    upstream_bp,
                    downstream_bp,
                    RepeatEnvironmentGeometryMode::CdsStopContext,
                    "annotated CDS stop-side context",
                    Some(context.transcript_id.clone()),
                )
            } else {
                Self::unavailable_repeat_geometry(
                    RepeatEnvironmentGeometryMode::CdsStopContext,
                    "no annotated CDS stop context available",
                )
            });
        } else {
            windows.push(Self::unavailable_repeat_geometry(
                RepeatEnvironmentGeometryMode::Transcript5utrStart,
                "no transcript projection available",
            ));
            windows.push(Self::unavailable_repeat_geometry(
                RepeatEnvironmentGeometryMode::Pol2PromoterUpstream,
                "no transcript projection available",
            ));
            windows.push(Self::unavailable_repeat_geometry(
                RepeatEnvironmentGeometryMode::CdsStopContext,
                "no transcript projection available",
            ));
        }
        windows
    }

    fn projected_repeat_transcript_contexts(
        repeat: &RepeatAnnotationRecord,
        transcripts: &[GenomeTranscriptRecord],
    ) -> Vec<RepeatTranscriptContext> {
        let context_span_start = repeat.start_1based.saturating_sub(10_000).max(1);
        let context_span_end = repeat.end_1based.saturating_add(10_000);
        let mut contexts = transcripts
            .iter()
            .filter(|transcript| {
                transcript
                    .chromosome
                    .eq_ignore_ascii_case(&repeat.chromosome)
            })
            .filter(|transcript| {
                transcript.transcript_start_1based <= context_span_end
                    && transcript.transcript_end_1based >= context_span_start
            })
            .map(|transcript| Self::transcript_context_for_repeat(repeat, transcript))
            .collect::<Vec<_>>();
        contexts.sort_by(|left, right| {
            let left_overlap_rank = if left.overlaps_repeat { 0 } else { 1 };
            let right_overlap_rank = if right.overlaps_repeat { 0 } else { 1 };
            left_overlap_rank
                .cmp(&right_overlap_rank)
                .then_with(|| {
                    left.signed_tss_distance_bp
                        .unwrap_or(i64::MAX)
                        .abs()
                        .cmp(&right.signed_tss_distance_bp.unwrap_or(i64::MAX).abs())
                })
                .then(left.transcript_id.cmp(&right.transcript_id))
        });
        contexts.truncate(12);
        contexts
    }

    fn load_all_genome_transcripts_for_repeat_projection(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<Vec<GenomeTranscriptRecord>, String> {
        let effective_catalog_path =
            catalog_path.unwrap_or(crate::genomes::default_catalog_discovery_token(false));
        let (catalog, _) = Self::open_reference_genome_catalog(Some(effective_catalog_path))
            .map_err(|e| e.to_string())?;
        let genes = catalog.list_gene_regions(genome_id, cache_dir)?;
        let mut transcripts = vec![];
        let mut seen = HashSet::new();
        for gene in genes {
            let gene_transcripts = catalog.list_gene_transcript_records(
                genome_id,
                &gene.chromosome,
                gene.start_1based,
                gene.end_1based,
                gene.gene_id.as_deref(),
                gene.gene_name.as_deref(),
                cache_dir,
            )?;
            for transcript in gene_transcripts {
                let key = format!("{}:{}", transcript.chromosome, transcript.transcript_id);
                if seen.insert(key) {
                    transcripts.push(transcript);
                }
            }
        }
        Ok(transcripts)
    }

    pub(crate) fn build_repeat_environment_cohort(
        &self,
        genome_id: &str,
        rmsk_path: &str,
        filter: &RepeatAnnotationFilter,
        upstream_bp: usize,
        downstream_bp: usize,
        geometry_mode: RepeatEnvironmentGeometryMode,
        limit: Option<usize>,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<RepeatEnvironmentCohortReport, EngineError> {
        let query = self.query_repeat_annotations(genome_id, rmsk_path, filter, limit)?;
        let mut warnings = query.warnings.clone();
        let transcripts = match Self::load_all_genome_transcripts_for_repeat_projection(
            genome_id,
            catalog_path,
            cache_dir,
        ) {
            Ok(transcripts) => transcripts,
            Err(err) => {
                warnings.push(format!(
                    "Transcript projection unavailable for repeat cohort '{}': {}",
                    genome_id, err
                ));
                vec![]
            }
        };
        let rows = query
            .rows
            .iter()
            .enumerate()
            .map(|(idx, repeat)| {
                let contexts = Self::projected_repeat_transcript_contexts(repeat, &transcripts);
                let geometry_windows = Self::build_repeat_geometry_windows(
                    repeat,
                    &contexts,
                    upstream_bp,
                    downstream_bp,
                );
                let selected_window = geometry_windows
                    .iter()
                    .find(|window| window.mode == geometry_mode)
                    .cloned()
                    .unwrap_or_else(|| {
                        Self::unavailable_repeat_geometry(
                            geometry_mode,
                            "requested geometry was not produced",
                        )
                    });
                let primary = Self::choose_primary_repeat_transcript_context(&contexts);
                let primary_gene_label = primary.map(|context| context.gene_label.clone());
                let primary_transcript_id = primary.map(|context| context.transcript_id.clone());
                RepeatEnvironmentCohortRow {
                    row_id: format!("repeat_{}", idx + 1),
                    repeat: repeat.clone(),
                    selected_geometry: geometry_mode,
                    selected_window,
                    geometry_windows,
                    transcript_contexts: contexts,
                    primary_gene_label,
                    primary_transcript_id,
                    rna_rank_score: None,
                    rna_support_labels: vec![],
                }
            })
            .collect::<Vec<_>>();
        Ok(RepeatEnvironmentCohortReport {
            schema: REPEAT_ENVIRONMENT_COHORT_REPORT_SCHEMA.to_string(),
            genome_id: genome_id.to_string(),
            rmsk_path: rmsk_path.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            filter: filter.clone(),
            selected_geometry: geometry_mode,
            upstream_bp,
            downstream_bp,
            parsed_repeat_count: query.parsed_row_count,
            matched_repeat_count: query.matched_row_count,
            returned_row_count: rows.len(),
            rows,
            warnings,
        })
    }

    fn reverse_complement_text(sequence: &str) -> String {
        sequence
            .as_bytes()
            .iter()
            .rev()
            .map(|base| IupacCode::letter_complement(*base) as char)
            .collect()
    }

    fn repeat_window_peak_genomic_position(
        window: &RepeatEnvironmentGeometryWindow,
        peak_position_0based: Option<usize>,
    ) -> Option<usize> {
        let peak = peak_position_0based?;
        let start = window.start_1based?;
        let end = window.end_1based?;
        if window.strand.as_deref() == Some("-") {
            Some(end.saturating_sub(peak))
        } else {
            Some(start.saturating_add(peak))
        }
    }

    fn repeat_tfbs_positive_support_fraction(track: &TfbsScoreTrackRow) -> f64 {
        if track.scored_window_count == 0 {
            return 0.0;
        }
        let positive_windows = track
            .forward_scores
            .iter()
            .zip(track.reverse_scores.iter())
            .filter(|(forward, reverse)| **forward > 0.0 || **reverse > 0.0)
            .count();
        positive_windows as f64 / track.scored_window_count as f64
    }

    pub(crate) fn summarize_window_cohort_tfbs(
        &self,
        cohort: &RepeatEnvironmentCohortReport,
        motifs: &[String],
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<WindowCohortTfbsReport, EngineError> {
        if motifs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SummarizeWindowCohortTfbs requires at least one motif".to_string(),
            });
        }
        let effective_catalog_path =
            catalog_path.unwrap_or(crate::genomes::default_catalog_discovery_token(false));
        let (catalog, _) = Self::open_reference_genome_catalog(Some(effective_catalog_path))?;
        let mut windows = vec![];
        let mut summary_rows = vec![];
        let mut warnings = cohort.warnings.clone();
        for row in &cohort.rows {
            let window = &row.selected_window;
            if !window.available {
                warnings.push(format!(
                    "Skipping repeat cohort row '{}' because {} geometry is unavailable: {}",
                    row.row_id,
                    cohort.selected_geometry.as_str(),
                    window.reason
                ));
                continue;
            }
            let (Some(chromosome), Some(start), Some(end)) = (
                window.chromosome.as_deref(),
                window.start_1based,
                window.end_1based,
            ) else {
                warnings.push(format!(
                    "Skipping repeat cohort row '{}' because selected window coordinates are incomplete",
                    row.row_id
                ));
                continue;
            };
            let sequence = match catalog.get_sequence_region_with_cache(
                &cohort.genome_id,
                chromosome,
                start,
                end,
                cache_dir,
            ) {
                Ok(sequence) => sequence,
                Err(err) => {
                    warnings.push(format!(
                        "Skipping repeat cohort row '{}' because {}:{}-{} could not be loaded: {}",
                        row.row_id, chromosome, start, end, err
                    ));
                    continue;
                }
            };
            let mut sequence_text = sequence;
            if window.strand.as_deref() == Some("-") {
                sequence_text = Self::reverse_complement_text(&sequence_text);
            }
            let label = row
                .primary_gene_label
                .as_deref()
                .unwrap_or(&row.repeat.normalized_alias)
                .to_string();
            let mut tfbs_score_tracks = self.summarize_tfbs_score_tracks(
                SequenceScanTarget::InlineSequence {
                    sequence_text,
                    topology: InlineSequenceTopology::Linear,
                    id_hint: Some(format!("{} {}", row.row_id, label)),
                    span_start_0based: None,
                    span_end_0based_exclusive: None,
                },
                motifs,
                score_kind,
                clip_negative,
            )?;
            if let Some(anchor) = window.anchor_1based {
                let position_0based = if window.strand.as_deref() == Some("-") {
                    end.saturating_sub(anchor)
                } else {
                    anchor.saturating_sub(start)
                };
                tfbs_score_tracks.tss_markers = vec![TfbsScoreTrackTssMarker {
                    feature_id: usize::MAX,
                    feature_kind: "repeat_environment_anchor".to_string(),
                    label: window.mode.as_str().to_string(),
                    position_0based,
                    is_reverse: window.strand.as_deref() == Some("-"),
                }];
            }
            for track in &tfbs_score_tracks.tracks {
                summary_rows.push(WindowCohortTfbsSummaryRow {
                    row_id: row.row_id.clone(),
                    label: label.clone(),
                    tf_id: track.tf_id.clone(),
                    tf_name: track.tf_name.clone(),
                    max_score: track.max_score,
                    peak_position_0based: track.max_position_0based,
                    peak_genomic_position_1based: Self::repeat_window_peak_genomic_position(
                        window,
                        track.max_position_0based,
                    ),
                    positive_fraction: Self::repeat_tfbs_positive_support_fraction(track),
                });
            }
            windows.push(WindowCohortTfbsWindowReport {
                row_id: row.row_id.clone(),
                label,
                chromosome: chromosome.to_string(),
                start_1based: start,
                end_1based: end,
                anchor_1based: window.anchor_1based,
                geometry_mode: window.mode,
                tfbs_score_tracks,
            });
        }
        Ok(WindowCohortTfbsReport {
            schema: WINDOW_COHORT_TFBS_REPORT_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            genome_id: cohort.genome_id.clone(),
            cohort_schema: cohort.schema.clone(),
            geometry_mode: cohort.selected_geometry,
            score_kind,
            clip_negative,
            motifs_requested: motifs.to_vec(),
            returned_window_count: windows.len(),
            windows,
            summary_rows,
            warnings,
        })
    }
}
