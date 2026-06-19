use super::*;

impl GentleEngine {
    pub(in crate::engine) fn interpret_probe_region_evidence(
        dna: &DNAsequence,
        seq_id: &str,
        gene_label: Option<&str>,
        level: Option<&str>,
        min_abs_logfc: Option<f64>,
    ) -> Result<ProbeRegionEvidenceInterpretationReport, EngineError> {
        if let Some(threshold) = min_abs_logfc
            && (!threshold.is_finite() || threshold < 0.0)
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "InterpretProbeRegionEvidence min_abs_logfc must be >= 0".to_string(),
                cause_chain: vec![],
            });
        }
        let level = Self::probe_region_interpretation_level(level)?;
        let gene_label = gene_label
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string);
        let mut warnings = Vec::new();
        let evidence = Self::projected_probe_region_evidence(dna, &level, min_abs_logfc);
        let transcripts = Self::probe_region_transcript_models(dna, gene_label.as_deref());

        if evidence.is_empty() {
            warnings.push(
                "No projected probe-region array features matched the requested interpretation filter"
                    .to_string(),
            );
        }
        if transcripts.is_empty() {
            warnings.push(
                "No transcript/exon models matched the requested gene filter; evidence is reported without transcript compatibility calls"
                    .to_string(),
            );
        }
        let (coordinate_frame, coordinate_system, coordinate_chromosome, coordinate_warnings) =
            Self::probe_region_evidence_report_coordinate_summary(&evidence);
        warnings.extend(coordinate_warnings);

        let mut transcript_counts = transcripts
            .iter()
            .map(|tx| {
                (
                    tx.transcript_id.clone(),
                    ProbeRegionTranscriptEvidenceCounts::default(),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let mut rows = Vec::new();
        for item in &evidence {
            let mut overlapping_transcripts = Vec::new();
            let mut overlapping_exon_count = 0usize;
            let mut span_overlaps = Vec::new();
            let mut transcript_mappings = Vec::new();
            for tx in &transcripts {
                if let Some(mapping) = Self::probe_region_evidence_transcript_mapping(item, tx) {
                    if !mapping.exon_ordinals.is_empty() {
                        overlapping_exon_count += mapping.exon_ordinals.len();
                        overlapping_transcripts.push(tx.transcript_id.clone());
                    } else {
                        span_overlaps.push(tx.transcript_id.clone());
                    }
                    transcript_mappings.push(mapping);
                }
            }
            overlapping_transcripts.sort();
            overlapping_transcripts.dedup();
            span_overlaps.sort();
            span_overlaps.dedup();
            transcript_mappings.sort_by(|left, right| {
                left.transcript_id
                    .cmp(&right.transcript_id)
                    .then(left.mapping_kind.cmp(&right.mapping_kind))
            });

            let mapping_status = if transcripts.is_empty() {
                "no_transcript_models"
            } else if overlapping_transcripts.len() > 1 {
                "shared_exon_overlap"
            } else if overlapping_transcripts.len() == 1 && overlapping_exon_count > 1 {
                "unique_multi_exon_overlap"
            } else if overlapping_transcripts.len() == 1 {
                "unique_exon_overlap"
            } else if !span_overlaps.is_empty() {
                "no_exon_overlap_inside_transcript_span"
            } else {
                "no_transcript_overlap"
            }
            .to_string();
            let relationship = if transcripts.is_empty() {
                "no_transcript_models"
            } else if !overlapping_transcripts.is_empty() {
                "compatible_with_exon_geometry"
            } else if !span_overlaps.is_empty() {
                "constrains_by_non_exonic_overlap"
            } else {
                "unmapped_to_transcript_models"
            }
            .to_string();

            let mut ambiguity_tags = BTreeSet::from([
                "multi_hit_not_assessed".to_string(),
                "probe_sequence_alignment_not_assessed".to_string(),
                "isoform_support_not_inferred".to_string(),
            ]);
            if overlapping_transcripts.len() > 1 {
                ambiguity_tags.insert("shared_transcript_overlap".to_string());
            }
            if item.parent_feature_id.is_some() {
                ambiguity_tags.insert("parent_probeset_context".to_string());
            }
            if item.intensity_source.as_deref() == Some("probe_level_input") {
                ambiguity_tags.insert("pm_probe_input".to_string());
            }
            if item.assembly_check.as_deref() == Some("projected_from_native_coordinate_system") {
                ambiguity_tags.insert("coordinate_projection_used".to_string());
            }
            if transcripts.is_empty() {
                ambiguity_tags.insert("no_transcript_models".to_string());
            }
            if !span_overlaps.is_empty() && overlapping_transcripts.is_empty() {
                ambiguity_tags.insert("non_exonic_transcript_span_overlap".to_string());
            }

            for mapping in &transcript_mappings {
                if let Some(counts) = transcript_counts.get_mut(&mapping.transcript_id) {
                    if mapping.exon_ordinals.is_empty() {
                        counts.constraining += 1;
                        counts.constraining_score += mapping.geometry_score;
                    } else {
                        counts.compatible += 1;
                        counts.compatible_score += mapping.geometry_score;
                        if overlapping_transcripts.len() > 1 {
                            counts.shared += 1;
                            counts.shared_score += mapping.geometry_score;
                        } else {
                            counts.unique += 1;
                            counts.unique_score += mapping.geometry_score;
                        }
                    }
                }
            }

            rows.push(ProbeRegionEvidenceMappingRow {
                evidence_id: item.evidence_id.clone(),
                level: item.level.clone(),
                feature_id: item.feature_id.clone(),
                parent_feature_id: item.parent_feature_id.clone(),
                intensity_source: item.intensity_source.clone(),
                chromosome: item.chromosome.clone(),
                start_1based: item.start_1based,
                end_1based: item.end_1based,
                strand: item.strand.clone(),
                logfc: item.logfc,
                overlapping_transcript_ids: overlapping_transcripts,
                overlapping_exon_count,
                transcript_mappings,
                mapping_status,
                ambiguity_tags: ambiguity_tags.into_iter().collect(),
                relationship,
            });
        }

        let total_evidence = evidence.len();
        let transcript_rows = transcripts
            .iter()
            .map(|tx| {
                let counts = transcript_counts
                    .remove(&tx.transcript_id)
                    .unwrap_or_default();
                let unmapped =
                    total_evidence.saturating_sub(counts.compatible + counts.constraining);
                let relationship_summary = if counts.unique > 0 {
                    "has_unique_compatible_evidence"
                } else if counts.shared > 0 {
                    "only_shared_compatible_evidence"
                } else if counts.constraining > 0 {
                    "has_non_exonic_constraining_evidence"
                } else {
                    "no_mapped_evidence"
                }
                .to_string();
                let review_status = if counts.unique > 0 {
                    "unique_geometry_for_review"
                } else if counts.shared > 0 {
                    "shared_geometry_for_review"
                } else if counts.constraining > 0 {
                    "non_exonic_constraint_for_review"
                } else {
                    "no_probe_region_evidence"
                }
                .to_string();
                ProbeRegionEvidenceTranscriptRow {
                    transcript_id: tx.transcript_id.clone(),
                    gene: tx.gene.clone(),
                    label: tx.label.clone(),
                    strand: tx.strand.clone(),
                    exon_count: tx.exon_ranges_0based.len(),
                    compatible_evidence_count: counts.compatible,
                    constraining_evidence_count: counts.constraining,
                    shared_evidence_count: counts.shared,
                    unique_evidence_count: counts.unique,
                    unmapped_evidence_count: unmapped,
                    compatible_geometry_score: counts.compatible_score,
                    shared_geometry_score: counts.shared_score,
                    unique_geometry_score: counts.unique_score,
                    constraining_geometry_score: counts.constraining_score,
                    review_status,
                    relationship_summary,
                }
            })
            .collect::<Vec<_>>();

        Ok(ProbeRegionEvidenceInterpretationReport {
            schema: PROBE_REGION_EVIDENCE_INTERPRETATION_SCHEMA.to_string(),
            seq_id: seq_id.to_string(),
            gene_label,
            level,
            coordinate_frame,
            coordinate_system,
            coordinate_chromosome,
            min_abs_logfc,
            array_feature_count: evidence.len(),
            transcript_count: transcripts.len(),
            evidence_rows: rows,
            transcript_rows,
            warnings,
        })
    }

    pub(super) fn probe_region_interpretation_level(
        level: Option<&str>,
    ) -> Result<String, EngineError> {
        let normalized = level
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("all")
            .replace('-', "_")
            .to_ascii_lowercase();
        match normalized.as_str() {
            "all" | "*" => Ok("all".to_string()),
            "probe_region" | "region" | "probeset" | "probeset_region" | "psr" => {
                Ok("probe_region".to_string())
            }
            "pm_probe" | "probe" | "probe_level" | "pm" => Ok("pm_probe".to_string()),
            other => Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "InterpretProbeRegionEvidence level '{other}' is not supported; use all, probe_region, or pm_probe"
                ),
                cause_chain: vec![],
            }),
        }
    }

    pub(super) fn projected_probe_region_evidence(
        dna: &DNAsequence,
        level: &str,
        min_abs_logfc: Option<f64>,
    ) -> Vec<ProbeRegionProjectedEvidence> {
        let mut out = Vec::new();
        for (idx, feature) in dna.features().iter().enumerate() {
            if Self::probe_region_first_qualifier(feature, &["gentle_track_source"]).as_deref()
                != Some("Array")
                || Self::probe_region_first_qualifier(feature, &["gentle_array_dataset"]).as_deref()
                    != Some("probe_region_output")
            {
                continue;
            }
            let feature_level =
                Self::probe_region_first_qualifier(feature, &["gentle_array_level"])
                    .unwrap_or_else(|| "probe_region".to_string());
            if level != "all" && feature_level != level {
                continue;
            }
            let logfc = Self::probe_region_first_qualifier(feature, &["logFC", "score"]).and_then(
                |value| {
                    value
                        .trim()
                        .parse::<f64>()
                        .ok()
                        .filter(|value| value.is_finite())
                },
            );
            if let Some(threshold) = min_abs_logfc
                && logfc.map(|value| value.abs() < threshold).unwrap_or(true)
            {
                continue;
            }
            let mut ranges = Vec::new();
            collect_location_ranges_usize(&feature.location, &mut ranges);
            if ranges.is_empty() {
                continue;
            }
            let feature_id = Self::probe_region_first_qualifier(
                feature,
                &["gentle_array_feature_id", "feature_id", "label"],
            )
            .unwrap_or_else(|| format!("array_feature_{idx}"));
            let contrast = Self::probe_region_first_qualifier(feature, &["gentle_array_contrast"]);
            let evidence_id = contrast
                .as_deref()
                .map(|contrast| format!("{feature_id}:{contrast}"))
                .unwrap_or_else(|| feature_id.clone());
            let (local_start, local_end) = Self::probe_region_range_span(&ranges).unwrap_or((0, 0));
            let genomic_start =
                Self::probe_region_first_qualifier(feature, &["genomic_start_1based"])
                    .and_then(|value| value.parse::<usize>().ok());
            let genomic_end = Self::probe_region_first_qualifier(feature, &["genomic_end_1based"])
                .and_then(|value| value.parse::<usize>().ok());
            let generic_start =
                Self::probe_region_first_qualifier(feature, &["start_1based", "start"])
                    .and_then(|value| value.parse::<usize>().ok());
            let generic_end =
                Self::probe_region_first_qualifier(feature, &["end_1based", "stop", "end"])
                    .and_then(|value| value.parse::<usize>().ok());
            let anchor_genome_id =
                Self::probe_region_first_qualifier(feature, &["gentle_array_anchor_genome_id"]);
            let anchor_start_1based =
                Self::probe_region_first_qualifier(feature, &["gentle_array_anchor_start_1based"])
                    .and_then(|value| value.parse::<usize>().ok());
            let anchor_end_1based =
                Self::probe_region_first_qualifier(feature, &["gentle_array_anchor_end_1based"])
                    .and_then(|value| value.parse::<usize>().ok());
            let anchor_strand =
                Self::probe_region_first_qualifier(feature, &["gentle_array_anchor_strand"])
                    .and_then(|value| value.chars().next())
                    .filter(|value| matches!(value, '+' | '-'));
            let coordinate_frame = if genomic_start.is_some()
                || genomic_end.is_some()
                || anchor_start_1based.is_some()
                || anchor_end_1based.is_some()
            {
                "genomic_1based"
            } else {
                "sequence_local_1based"
            }
            .to_string();
            out.push(ProbeRegionProjectedEvidence {
                evidence_id,
                level: feature_level,
                feature_id,
                parent_feature_id: Self::probe_region_first_qualifier(
                    feature,
                    &["gentle_array_parent_feature_id"],
                ),
                intensity_source: Self::probe_region_first_qualifier(
                    feature,
                    &["gentle_array_intensity_source"],
                ),
                chromosome: Self::probe_region_first_qualifier(feature, &["chromosome"]),
                start_1based: genomic_start
                    .or(generic_start)
                    .or_else(|| Some(local_start + 1)),
                end_1based: genomic_end.or(generic_end).or_else(|| Some(local_end)),
                strand: Self::probe_region_first_qualifier(feature, &["strand", "array_strand"]),
                logfc,
                assembly_check: Self::probe_region_first_qualifier(
                    feature,
                    &["gentle_array_assembly_check"],
                ),
                coordinate_frame,
                anchor_genome_id,
                anchor_start_1based,
                anchor_end_1based,
                anchor_strand,
                ranges_0based: ranges,
            });
        }
        out
    }

    pub(super) fn probe_region_transcript_models(
        dna: &DNAsequence,
        gene_label: Option<&str>,
    ) -> Vec<ProbeRegionTranscriptModel> {
        let mut exon_ranges_by_transcript: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
        for feature in dna.features() {
            if !feature.kind.eq_ignore_ascii_case("exon") {
                continue;
            }
            let mut ranges = Vec::new();
            collect_location_ranges_usize(&feature.location, &mut ranges);
            if ranges.is_empty() {
                continue;
            }
            for transcript_id in Self::probe_region_feature_text_values(
                feature,
                &["transcript_id", "Parent", "parent", "transcript"],
            ) {
                exon_ranges_by_transcript
                    .entry(transcript_id)
                    .or_default()
                    .extend(ranges.iter().copied());
            }
        }

        let mut out = Vec::new();
        let mut seen = BTreeSet::new();
        for (idx, feature) in dna.features().iter().enumerate() {
            if !Self::probe_region_is_transcript_feature(feature) {
                continue;
            }
            if !Self::probe_region_feature_matches_gene(feature, gene_label) {
                continue;
            }
            let transcript_id = Self::probe_region_first_qualifier(
                feature,
                &[
                    "transcript_id",
                    "transcript",
                    "transcript_cluster_id",
                    "ID",
                    "id",
                    "Name",
                    "label",
                ],
            )
            .unwrap_or_else(|| format!("transcript_{idx}"));
            if !seen.insert(transcript_id.clone()) {
                continue;
            }
            let mut ranges = Vec::new();
            collect_location_ranges_usize(&feature.location, &mut ranges);
            if ranges.len() <= 1
                && let Some(exon_ranges) = exon_ranges_by_transcript.get(&transcript_id)
            {
                ranges = exon_ranges.clone();
            }
            if ranges.is_empty() {
                continue;
            }
            ranges.sort();
            ranges.dedup();
            let strand = Self::probe_region_first_qualifier(feature, &["strand"]);
            let is_reverse = strand.as_deref() == Some("-");
            let exon_count = ranges.len();
            let exons = ranges
                .iter()
                .enumerate()
                .map(|(idx, range)| ProbeRegionTranscriptExon {
                    ordinal: if is_reverse {
                        exon_count.saturating_sub(idx)
                    } else {
                        idx.saturating_add(1)
                    },
                    range_0based: *range,
                })
                .collect::<Vec<_>>();
            out.push(ProbeRegionTranscriptModel {
                transcript_id,
                gene: Self::probe_region_first_qualifier(
                    feature,
                    &["gene", "gene_symbol", "gene_name", "gene_id"],
                ),
                label: Self::probe_region_first_qualifier(feature, &["label", "Name"]),
                strand,
                exons,
                span_0based: Self::probe_region_range_span(&ranges),
                exon_ranges_0based: ranges,
            });
        }
        out
    }

    pub(super) fn probe_region_is_transcript_feature(feature: &gb_io::seq::Feature) -> bool {
        matches!(
            feature.kind.to_ascii_lowercase().as_str(),
            "mrna" | "transcript" | "ncrna" | "misc_rna" | "rrna" | "trna"
        )
    }

    pub(super) fn probe_region_feature_matches_gene(
        feature: &gb_io::seq::Feature,
        gene_label: Option<&str>,
    ) -> bool {
        let Some(requested) = gene_label.map(str::trim).filter(|value| !value.is_empty()) else {
            return true;
        };
        let requested_lower = requested.to_ascii_lowercase();
        let values = Self::probe_region_feature_text_values(
            feature,
            &[
                "gene",
                "gene_symbol",
                "gene_name",
                "gene_id",
                "label",
                "Name",
                "transcript_id",
                "transcript",
            ],
        );
        values.iter().any(|value| {
            let lower = value.to_ascii_lowercase();
            lower == requested_lower || lower.contains(&requested_lower)
        })
    }

    pub(super) fn probe_region_feature_text_values(
        feature: &gb_io::seq::Feature,
        keys: &[&str],
    ) -> Vec<String> {
        let mut out = Vec::new();
        let mut seen = BTreeSet::new();
        for key in keys {
            for value in feature.qualifier_values(key) {
                for part in value
                    .split([',', ';'])
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    if seen.insert(part.to_ascii_lowercase()) {
                        out.push(part.to_string());
                    }
                }
            }
        }
        out
    }

    pub(super) fn probe_region_first_qualifier(
        feature: &gb_io::seq::Feature,
        keys: &[&str],
    ) -> Option<String> {
        keys.iter().find_map(|key| {
            feature
                .qualifier_values(key)
                .find(|value| !value.trim().is_empty())
                .map(|value| value.trim().to_string())
        })
    }

    pub(super) fn probe_region_ranges_overlap(
        left: &(usize, usize),
        right: &(usize, usize),
    ) -> bool {
        left.0 < right.1 && right.0 < left.1
    }

    pub(super) fn probe_region_range_overlap_bp(
        left: &(usize, usize),
        right: &(usize, usize),
    ) -> usize {
        let start = left.0.max(right.0);
        let end = left.1.min(right.1);
        end.saturating_sub(start)
    }

    pub(super) fn probe_region_range_span(ranges: &[(usize, usize)]) -> Option<(usize, usize)> {
        let start = ranges.iter().map(|(start, _)| *start).min()?;
        let end = ranges.iter().map(|(_, end)| *end).max()?;
        Some((start, end))
    }

    pub(super) fn probe_region_local_range_string(range: &(usize, usize)) -> String {
        format!("{}..{}", range.0.saturating_add(1), range.1)
    }

    pub(super) fn probe_region_evidence_report_coordinate_summary(
        evidence: &[ProbeRegionProjectedEvidence],
    ) -> (String, Option<String>, Option<String>, Vec<String>) {
        let frames = evidence
            .iter()
            .map(|item| item.coordinate_frame.clone())
            .collect::<BTreeSet<_>>();
        let coordinate_frame = match frames.len() {
            0 => "unknown".to_string(),
            1 => frames
                .into_iter()
                .next()
                .unwrap_or_else(|| "unknown".to_string()),
            _ => "mixed".to_string(),
        };
        let systems = evidence
            .iter()
            .filter_map(|item| item.anchor_genome_id.clone())
            .collect::<BTreeSet<_>>();
        let coordinate_system = (systems.len() == 1).then(|| {
            systems
                .iter()
                .next()
                .cloned()
                .unwrap_or_else(|| "unknown".to_string())
        });
        let chromosomes = evidence
            .iter()
            .filter_map(|item| item.chromosome.clone())
            .collect::<BTreeSet<_>>();
        let coordinate_chromosome = (chromosomes.len() == 1).then(|| {
            chromosomes
                .iter()
                .next()
                .cloned()
                .unwrap_or_else(|| "unknown".to_string())
        });
        let mut warnings = Vec::new();
        if coordinate_frame == "mixed" {
            warnings.push(
                "probe_region_evidence_report_mixes_coordinate_frames; renderer uses row-level coordinates without local realignment".to_string(),
            );
        }
        if systems.len() > 1 {
            warnings.push(
                "probe_region_evidence_report_mixes_coordinate_systems; review coordinates before publication use".to_string(),
            );
        }
        if chromosomes.len() > 1 {
            warnings.push(
                "probe_region_evidence_report_mixes_chromosomes; review coordinates before publication use".to_string(),
            );
        }
        (
            coordinate_frame,
            coordinate_system,
            coordinate_chromosome,
            warnings,
        )
    }

    pub(super) fn probe_region_evidence_local_interval_to_report_span(
        item: &ProbeRegionProjectedEvidence,
        range: &(usize, usize),
    ) -> (usize, usize) {
        if item.coordinate_frame == "genomic_1based" {
            if let (Some(anchor_start), Some(anchor_end)) =
                (item.anchor_start_1based, item.anchor_end_1based)
            {
                if item.anchor_strand == Some('-') {
                    let start = anchor_end.saturating_sub(range.1.saturating_sub(1));
                    let end = anchor_end.saturating_sub(range.0);
                    return (start.min(end).max(1), start.max(end).max(1));
                }
                let start = anchor_start.saturating_add(range.0);
                let end = anchor_start.saturating_add(range.1.saturating_sub(1));
                return (start.min(end).max(1), start.max(end).max(1));
            }
            if let (Some((local_start, _)), Some(evidence_start)) = (
                Self::probe_region_range_span(&item.ranges_0based),
                item.start_1based,
            ) {
                let report_start = evidence_start as i128 + range.0 as i128 - local_start as i128;
                let report_end = evidence_start as i128 + range.1.saturating_sub(1) as i128
                    - local_start as i128;
                let start = report_start.max(1) as usize;
                let end = report_end.max(1) as usize;
                return (start.min(end), start.max(end));
            }
        }
        (range.0.saturating_add(1), range.1)
    }

    pub(super) fn probe_region_evidence_local_interval_to_report_range_string(
        item: &ProbeRegionProjectedEvidence,
        range: &(usize, usize),
    ) -> String {
        let (start, end) = Self::probe_region_evidence_local_interval_to_report_span(item, range);
        format!("{start}..{end}")
    }

    pub(super) fn probe_region_evidence_transcript_mapping(
        item: &ProbeRegionProjectedEvidence,
        tx: &ProbeRegionTranscriptModel,
    ) -> Option<ProbeRegionEvidenceTranscriptMapping> {
        let mut exon_hits = tx
            .exons
            .iter()
            .filter_map(|exon| {
                let overlap_bp = item
                    .ranges_0based
                    .iter()
                    .map(|evidence_range| {
                        Self::probe_region_range_overlap_bp(&exon.range_0based, evidence_range)
                    })
                    .sum::<usize>();
                (overlap_bp > 0).then_some((exon, overlap_bp))
            })
            .collect::<Vec<_>>();
        exon_hits.sort_by(|left, right| left.0.ordinal.cmp(&right.0.ordinal));

        if !exon_hits.is_empty() {
            let exon_ordinals = exon_hits
                .iter()
                .map(|(exon, _)| exon.ordinal)
                .collect::<Vec<_>>();
            let local_exon_ranges_1based = exon_hits
                .iter()
                .map(|(exon, _)| Self::probe_region_local_range_string(&exon.range_0based))
                .collect::<Vec<_>>();
            let exon_ranges_1based = exon_hits
                .iter()
                .map(|(exon, _)| {
                    Self::probe_region_evidence_local_interval_to_report_range_string(
                        item,
                        &exon.range_0based,
                    )
                })
                .collect::<Vec<_>>();
            let overlap_bp = exon_hits.iter().map(|(_, overlap)| *overlap).sum();
            let junction_spans = Self::probe_region_evidence_junction_spans(item, tx, &exon_hits);
            let mapping_kind = if !junction_spans.is_empty() {
                "junction_spanning_exon_overlap"
            } else if exon_ordinals.len() > 1 {
                "multi_exon_overlap"
            } else {
                "exon_overlap"
            }
            .to_string();
            let (geometry_score, geometry_score_class, score_basis) =
                Self::probe_region_evidence_geometry_score(
                    &mapping_kind,
                    overlap_bp,
                    junction_spans.len(),
                );
            return Some(ProbeRegionEvidenceTranscriptMapping {
                transcript_id: tx.transcript_id.clone(),
                coordinate_frame: item.coordinate_frame.clone(),
                mapping_kind,
                geometry_score,
                geometry_score_class,
                score_basis,
                exon_ordinals,
                exon_ranges_1based,
                local_exon_ranges_1based,
                junction_spans,
                overlap_bp,
            });
        }

        let span = tx.span_0based?;
        let overlap_bp = item
            .ranges_0based
            .iter()
            .map(|range| Self::probe_region_range_overlap_bp(&span, range))
            .sum::<usize>();
        (overlap_bp > 0).then(|| ProbeRegionEvidenceTranscriptMapping {
            transcript_id: tx.transcript_id.clone(),
            coordinate_frame: item.coordinate_frame.clone(),
            mapping_kind: "transcript_span_non_exonic".to_string(),
            geometry_score: -0.25,
            geometry_score_class: "constraining_non_exonic_geometry".to_string(),
            score_basis: vec![
                "mapping_kind=transcript_span_non_exonic".to_string(),
                format!("overlap_bp={overlap_bp}"),
                "probe_sequence_alignment_not_assessed".to_string(),
                "isoform_support_not_inferred".to_string(),
            ],
            exon_ordinals: Vec::new(),
            exon_ranges_1based: Vec::new(),
            local_exon_ranges_1based: Vec::new(),
            junction_spans: Vec::new(),
            overlap_bp,
        })
    }

    pub(super) fn probe_region_evidence_geometry_score(
        mapping_kind: &str,
        overlap_bp: usize,
        junction_span_count: usize,
    ) -> (f64, String, Vec<String>) {
        let (score, class) = match mapping_kind {
            "junction_spanning_exon_overlap" => (0.75, "junction_spanning_geometry"),
            "multi_exon_overlap" => (0.60, "multi_exon_geometry"),
            "exon_overlap" => (0.50, "exon_geometry"),
            _ => (0.0, "unscored_geometry"),
        };
        (
            score,
            class.to_string(),
            vec![
                format!("mapping_kind={mapping_kind}"),
                format!("overlap_bp={overlap_bp}"),
                format!("junction_spans={junction_span_count}"),
                "probe_sequence_alignment_not_assessed".to_string(),
                "isoform_support_not_inferred".to_string(),
            ],
        )
    }

    pub(super) fn probe_region_evidence_junction_spans(
        item: &ProbeRegionProjectedEvidence,
        tx: &ProbeRegionTranscriptModel,
        exon_hits: &[(&ProbeRegionTranscriptExon, usize)],
    ) -> Vec<ProbeRegionEvidenceJunctionSpan> {
        let hit_ordinals = exon_hits
            .iter()
            .map(|(exon, _)| exon.ordinal)
            .collect::<BTreeSet<_>>();
        let mut out = Vec::new();
        let mut genomic_exons = tx.exons.iter().collect::<Vec<_>>();
        genomic_exons.sort_by(|left, right| {
            left.range_0based
                .0
                .cmp(&right.range_0based.0)
                .then(left.range_0based.1.cmp(&right.range_0based.1))
        });
        for pair in genomic_exons.windows(2) {
            let left = pair[0];
            let right = pair[1];
            if !hit_ordinals.contains(&left.ordinal) || !hit_ordinals.contains(&right.ordinal) {
                continue;
            }
            let gap = (left.range_0based.1, right.range_0based.0);
            let one_range_spans_boundary = item
                .ranges_0based
                .iter()
                .any(|range| range.0 < left.range_0based.1 && range.1 > right.range_0based.0);
            let joined_evidence_hits_both_sides = item.ranges_0based.len() > 1
                && item
                    .ranges_0based
                    .iter()
                    .any(|range| Self::probe_region_ranges_overlap(&left.range_0based, range))
                && item
                    .ranges_0based
                    .iter()
                    .any(|range| Self::probe_region_ranges_overlap(&right.range_0based, range));
            if !one_range_spans_boundary && !joined_evidence_hits_both_sides {
                continue;
            }
            let is_reverse = tx.strand.as_deref() == Some("-");
            let (from_exon_ordinal, to_exon_ordinal) = if is_reverse {
                (right.ordinal, left.ordinal)
            } else {
                (left.ordinal, right.ordinal)
            };
            let (report_gap_start, report_gap_end) =
                Self::probe_region_evidence_local_interval_to_report_span(item, &gap);
            out.push(ProbeRegionEvidenceJunctionSpan {
                from_exon_ordinal,
                to_exon_ordinal,
                genomic_start_1based: report_gap_start,
                genomic_end_1based: report_gap_end,
                local_start_1based: Some(gap.0.saturating_add(1)),
                local_end_1based: Some(gap.1),
            });
        }
        out.sort_by(|left, right| {
            left.from_exon_ordinal
                .cmp(&right.from_exon_ordinal)
                .then(left.to_exon_ordinal.cmp(&right.to_exon_ordinal))
        });
        out
    }
}
