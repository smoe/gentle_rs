//! Shared promoter-window, promoter-context, and allele-materialization helpers.
//!
//! These routines back the ClawBio-facing promoter-SNP-to-luciferase flow and
//! are kept in the engine so GUI/CLI/shell/workflow paths can reuse the same
//! deterministic behavior.

use super::*;

impl GentleEngine {
    fn feature_identifier_values(
        feature: &gb_io::seq::Feature,
        feature_id: usize,
        extra_keys: &[&str],
    ) -> Vec<String> {
        let mut out = vec![
            Self::feature_display_label(feature, feature_id),
            format!("{}", feature_id + 1),
            format!("#{}", feature_id + 1),
        ];
        let mut seen = out
            .iter()
            .map(|value| value.to_ascii_uppercase())
            .collect::<HashSet<_>>();
        let mut push_value = |raw: &str| {
            let trimmed = raw.trim();
            if trimmed.is_empty() {
                return;
            }
            let normalized = trimmed.to_ascii_uppercase();
            if seen.insert(normalized) {
                out.push(trimmed.to_string());
            }
        };
        for key in [
            "label",
            "name",
            "standard_name",
            "gene",
            "gene_id",
            "gene_name",
            "transcript_id",
            "locus_tag",
            "product",
            "note",
            "db_xref",
        ] {
            for value in feature.qualifier_values(key) {
                push_value(value);
            }
        }
        for key in extra_keys {
            for value in feature.qualifier_values(key) {
                push_value(value);
            }
        }
        out
    }

    fn feature_matches_identifier(
        feature: &gb_io::seq::Feature,
        feature_id: usize,
        query: &str,
        extra_keys: &[&str],
    ) -> bool {
        let trimmed = query.trim();
        if trimmed.is_empty() {
            return false;
        }
        Self::feature_identifier_values(feature, feature_id, extra_keys)
            .into_iter()
            .any(|candidate| {
                candidate.eq_ignore_ascii_case(trimmed)
                    || candidate
                        .rsplit(':')
                        .next()
                        .is_some_and(|tail| tail.eq_ignore_ascii_case(trimmed))
            })
    }

    fn feature_span_bounds(feature: &gb_io::seq::Feature) -> Option<(usize, usize)> {
        let mut ranges = vec![];
        collect_location_ranges_usize(&feature.location, &mut ranges);
        if ranges.is_empty() {
            None
        } else {
            let start = ranges.iter().map(|(start, _)| *start).min()?;
            let end = ranges.iter().map(|(_, end)| *end).max()?;
            (end > start).then_some((start, end))
        }
    }

    fn transcript_gene_metadata(
        dna: &DNAsequence,
        transcript_feature: &gb_io::seq::Feature,
        _transcript_feature_id: usize,
    ) -> (Option<String>, Option<String>) {
        let direct_gene_label = Self::first_nonempty_feature_qualifier(
            transcript_feature,
            &["gene", "gene_name", "locus_tag"],
        );
        let direct_gene_id =
            Self::first_nonempty_feature_qualifier(transcript_feature, &["gene_id", "locus_tag"]);
        if direct_gene_label.is_some() || direct_gene_id.is_some() {
            return (direct_gene_label, direct_gene_id);
        }

        let Some((tx_start, tx_end)) = Self::feature_span_bounds(transcript_feature) else {
            return (None, None);
        };
        let tx_reverse = feature_is_reverse(transcript_feature);
        let mut candidates = dna
            .features()
            .iter()
            .enumerate()
            .filter_map(|(feature_id, feature)| {
                (Self::construct_reasoning_role_from_feature(feature) == Some(ConstructRole::Gene)
                    && feature_is_reverse(feature) == tx_reverse
                    && Self::feature_overlaps_span(feature, tx_start, tx_end))
                .then(|| {
                    (
                        Self::feature_display_label(feature, feature_id),
                        Self::first_nonempty_feature_qualifier(feature, &["gene_id", "locus_tag"]),
                    )
                })
            })
            .collect::<Vec<_>>();
        candidates.sort_by(|left, right| {
            left.0
                .to_ascii_lowercase()
                .cmp(&right.0.to_ascii_lowercase())
                .then_with(|| left.1.cmp(&right.1))
        });
        candidates.into_iter().next().map_or((None, None), |row| {
            let gene_label = (!row.0.trim().is_empty()).then_some(row.0);
            (gene_label, row.1)
        })
    }

    fn collapse_promoter_window_records_by_gene(
        records: Vec<PromoterWindowRecord>,
    ) -> Vec<PromoterWindowRecord> {
        let mut grouped: BTreeMap<(String, String), Vec<PromoterWindowRecord>> = BTreeMap::new();
        for record in records {
            let key = (
                record
                    .gene_label
                    .clone()
                    .unwrap_or_else(|| record.transcript_id.clone())
                    .to_ascii_lowercase(),
                record.strand.clone(),
            );
            grouped.entry(key).or_default().push(record);
        }
        let mut collapsed = grouped
            .into_values()
            .filter_map(|mut group| {
                group.sort_by(|left, right| {
                    left.transcript_id
                        .to_ascii_lowercase()
                        .cmp(&right.transcript_id.to_ascii_lowercase())
                        .then_with(|| left.start_0based.cmp(&right.start_0based))
                        .then_with(|| left.end_0based_exclusive.cmp(&right.end_0based_exclusive))
                });
                let union_start = group.iter().map(|row| row.start_0based).min()?;
                let union_end = group.iter().map(|row| row.end_0based_exclusive).max()?;
                let mut representative = group.into_iter().next()?;
                representative.start_0based = union_start;
                representative.end_0based_exclusive = union_end;
                representative.source = "transcript_tss_gene_collapse".to_string();
                Some(representative)
            })
            .collect::<Vec<_>>();
        collapsed.sort_by(|left, right| {
            left.gene_label
                .clone()
                .unwrap_or_else(|| left.transcript_label.clone())
                .to_ascii_lowercase()
                .cmp(
                    &right
                        .gene_label
                        .clone()
                        .unwrap_or_else(|| right.transcript_label.clone())
                        .to_ascii_lowercase(),
                )
                .then_with(|| left.start_0based.cmp(&right.start_0based))
        });
        collapsed
    }

    fn collapse_promoter_window_records_by_exact_span(
        records: Vec<PromoterWindowRecord>,
    ) -> Vec<PromoterWindowRecord> {
        let mut grouped: BTreeMap<
            (
                Option<String>,
                Option<String>,
                String,
                usize,
                usize,
                usize,
                usize,
            ),
            Vec<PromoterWindowRecord>,
        > = BTreeMap::new();
        for record in records {
            let key = (
                record
                    .gene_label
                    .clone()
                    .map(|value| value.to_ascii_lowercase()),
                record
                    .gene_id
                    .clone()
                    .map(|value| value.to_ascii_lowercase()),
                record.strand.clone(),
                record.start_0based,
                record.end_0based_exclusive,
                record.upstream_bp,
                record.downstream_bp,
            );
            grouped.entry(key).or_default().push(record);
        }
        let mut collapsed = grouped
            .into_values()
            .filter_map(|mut group| {
                group.sort_by(|left, right| {
                    left.transcript_id
                        .to_ascii_lowercase()
                        .cmp(&right.transcript_id.to_ascii_lowercase())
                        .then_with(|| {
                            left.transcript_label
                                .to_ascii_lowercase()
                                .cmp(&right.transcript_label.to_ascii_lowercase())
                        })
                });
                let mut group_iter = group.into_iter();
                let mut representative = group_iter.next()?;
                let mut transcript_ids = representative.transcript_ids.clone();
                let mut transcript_labels = representative.transcript_labels.clone();
                let mut transcript_feature_id = representative.transcript_feature_id;
                for record in group_iter {
                    if transcript_feature_id.is_none() {
                        transcript_feature_id = record.transcript_feature_id;
                    }
                    transcript_ids.extend(record.transcript_ids);
                    transcript_labels.extend(record.transcript_labels);
                }
                transcript_ids.sort();
                transcript_ids.dedup_by(|left, right| left.eq_ignore_ascii_case(right));
                transcript_labels.sort();
                transcript_labels.dedup_by(|left, right| left.eq_ignore_ascii_case(right));
                representative.transcript_feature_id = transcript_feature_id;
                representative.transcript_count = transcript_ids.len().max(1);
                representative.transcript_ids = transcript_ids;
                representative.transcript_labels = transcript_labels;
                representative.source = "transcript_tss_exact_collapse".to_string();
                Some(representative)
            })
            .collect::<Vec<_>>();
        collapsed.sort_by(|left, right| {
            left.gene_label
                .clone()
                .unwrap_or_else(|| left.transcript_label.clone())
                .to_ascii_lowercase()
                .cmp(
                    &right
                        .gene_label
                        .clone()
                        .unwrap_or_else(|| right.transcript_label.clone())
                        .to_ascii_lowercase(),
                )
                .then_with(|| left.start_0based.cmp(&right.start_0based))
                .then_with(|| left.end_0based_exclusive.cmp(&right.end_0based_exclusive))
                .then_with(|| left.strand.cmp(&right.strand))
        });
        collapsed
    }

    pub(crate) fn derive_promoter_window_records(
        &self,
        dna: &DNAsequence,
        gene_label: Option<&str>,
        transcript_id: Option<&str>,
        upstream_bp: usize,
        downstream_bp: usize,
        collapse_mode: PromoterWindowCollapseMode,
    ) -> Vec<PromoterWindowRecord> {
        let seq_len = dna.len();
        let gene_label = gene_label.map(str::trim).filter(|value| !value.is_empty());
        let transcript_id = transcript_id
            .map(str::trim)
            .filter(|value| !value.is_empty());
        let mut records = dna
            .features()
            .iter()
            .enumerate()
            .filter_map(|(feature_id, feature)| {
                (Self::construct_reasoning_role_from_feature(feature)
                    == Some(ConstructRole::Transcript))
                .then_some((feature_id, feature))
            })
            .filter_map(|(feature_id, feature)| {
                if transcript_id.is_some_and(|query| {
                    !Self::feature_matches_identifier(
                        feature,
                        feature_id,
                        query,
                        &["transcript_id"],
                    )
                }) {
                    return None;
                }

                let (resolved_gene_label, resolved_gene_id) =
                    Self::transcript_gene_metadata(dna, feature, feature_id);
                if gene_label.is_some_and(|query| {
                    !resolved_gene_label
                        .as_deref()
                        .is_some_and(|value| value.eq_ignore_ascii_case(query))
                        && !Self::feature_matches_identifier(
                            feature,
                            feature_id,
                            query,
                            &["gene", "gene_id", "gene_name", "locus_tag"],
                        )
                }) {
                    return None;
                }

                let mut ranges = vec![];
                collect_location_ranges_usize(&feature.location, &mut ranges);
                if ranges.is_empty() {
                    return None;
                }
                ranges.sort_unstable_by(|left, right| {
                    left.0.cmp(&right.0).then(left.1.cmp(&right.1))
                });
                let is_reverse = feature_is_reverse(feature);
                let tss_local_0based = if is_reverse {
                    ranges
                        .iter()
                        .map(|(_, end)| end.saturating_sub(1))
                        .max()
                        .unwrap_or(0)
                } else {
                    ranges.iter().map(|(start, _)| *start).min().unwrap_or(0)
                };
                let (start_0based, end_0based_exclusive) = if is_reverse {
                    (
                        tss_local_0based.saturating_sub(downstream_bp),
                        tss_local_0based
                            .saturating_add(upstream_bp)
                            .saturating_add(1)
                            .min(seq_len),
                    )
                } else {
                    (
                        tss_local_0based.saturating_sub(upstream_bp),
                        tss_local_0based
                            .saturating_add(downstream_bp)
                            .saturating_add(1)
                            .min(seq_len),
                    )
                };
                if end_0based_exclusive <= start_0based {
                    return None;
                }
                let transcript_label = Self::feature_display_label(feature, feature_id);
                let transcript_id = Self::first_nonempty_feature_qualifier(
                    feature,
                    &["transcript_id", "name", "label", "product"],
                )
                .unwrap_or_else(|| transcript_label.clone());
                Some(PromoterWindowRecord {
                    gene_label: resolved_gene_label,
                    gene_id: resolved_gene_id,
                    transcript_id: transcript_id.clone(),
                    transcript_label: transcript_label.clone(),
                    transcript_feature_id: Some(feature_id),
                    transcript_count: 1,
                    transcript_ids: vec![transcript_id],
                    transcript_labels: vec![transcript_label],
                    strand: if is_reverse { "-" } else { "+" }.to_string(),
                    tss_local_0based,
                    start_0based,
                    end_0based_exclusive,
                    upstream_bp,
                    downstream_bp,
                    source: "transcript_tss".to_string(),
                })
            })
            .collect::<Vec<_>>();
        records.sort_by(|left, right| {
            left.gene_label
                .clone()
                .unwrap_or_else(|| left.transcript_label.clone())
                .to_ascii_lowercase()
                .cmp(
                    &right
                        .gene_label
                        .clone()
                        .unwrap_or_else(|| right.transcript_label.clone())
                        .to_ascii_lowercase(),
                )
                .then_with(|| {
                    left.transcript_id
                        .to_ascii_lowercase()
                        .cmp(&right.transcript_id.to_ascii_lowercase())
                })
                .then_with(|| left.start_0based.cmp(&right.start_0based))
                .then_with(|| left.end_0based_exclusive.cmp(&right.end_0based_exclusive))
        });
        records.dedup_by(|left, right| {
            left.transcript_id
                .eq_ignore_ascii_case(&right.transcript_id)
                && left.start_0based == right.start_0based
                && left.end_0based_exclusive == right.end_0based_exclusive
                && left.strand == right.strand
        });
        match collapse_mode {
            PromoterWindowCollapseMode::Transcript => records,
            PromoterWindowCollapseMode::Gene => {
                Self::collapse_promoter_window_records_by_gene(records)
            }
        }
    }

    fn build_promoter_window_feature(record: &PromoterWindowRecord) -> gb_io::seq::Feature {
        let base_location = gb_io::seq::Location::simple_range(
            record.start_0based as i64,
            record.end_0based_exclusive as i64,
        );
        let location = if record.strand == "-" {
            gb_io::seq::Location::Complement(Box::new(base_location))
        } else {
            base_location
        };
        let base_label = record
            .gene_label
            .clone()
            .unwrap_or_else(|| record.transcript_label.clone());
        let transcript_count = record
            .transcript_count
            .max(record.transcript_ids.len().max(1));
        let label = if transcript_count > 1 {
            format!("{base_label} promoter window ({transcript_count} tx)")
        } else {
            format!("{base_label} promoter window")
        };
        let note = if transcript_count > 1 {
            format!(
                "Generated promoter window from transcript TSS geometry shared by {transcript_count} transcripts (transcript_ids='{}', transcript_labels='{}', representative_tss_local_1based={}, upstream_bp={}, downstream_bp={}).",
                record.transcript_ids.join(","),
                record.transcript_labels.join(","),
                record.tss_local_0based.saturating_add(1),
                record.upstream_bp,
                record.downstream_bp
            )
        } else {
            format!(
                "Generated promoter window from transcript TSS geometry (transcript='{}', tss_local_1based={}, upstream_bp={}, downstream_bp={}).",
                record.transcript_label,
                record.tss_local_0based.saturating_add(1),
                record.upstream_bp,
                record.downstream_bp
            )
        };
        let mut qualifiers = vec![
            ("label".into(), Some(label)),
            ("note".into(), Some(note)),
            (
                "gentle_generated".into(),
                Some(ANNOTATE_PROMOTER_WINDOWS_GENERATED_TAG.to_string()),
            ),
            (
                "generated_by".into(),
                Some("AnnotatePromoterWindows".to_string()),
            ),
            ("promoter_source".into(), Some(record.source.clone())),
            ("regulatory_class".into(), Some("promoter".to_string())),
            ("transcript_id".into(), Some(record.transcript_id.clone())),
            (
                "transcript_count".into(),
                Some(transcript_count.to_string()),
            ),
            (
                "transcript_ids".into(),
                Some(record.transcript_ids.join(",")),
            ),
            (
                "transcript_label".into(),
                Some(record.transcript_label.clone()),
            ),
            (
                "transcript_labels".into(),
                Some(record.transcript_labels.join(",")),
            ),
            ("upstream_bp".into(), Some(record.upstream_bp.to_string())),
            (
                "downstream_bp".into(),
                Some(record.downstream_bp.to_string()),
            ),
        ];
        if let Some(gene_label) = record.gene_label.as_ref() {
            qualifiers.push(("gene".into(), Some(gene_label.clone())));
        }
        if let Some(gene_id) = record.gene_id.as_ref() {
            qualifiers.push(("gene_id".into(), Some(gene_id.clone())));
        }
        gb_io::seq::Feature {
            kind: "promoter".into(),
            location,
            qualifiers,
        }
    }

    pub(crate) fn annotate_promoter_windows_for_sequence(
        &mut self,
        seq_id: &str,
        gene_label: Option<&str>,
        transcript_id: Option<&str>,
        upstream_bp: usize,
        downstream_bp: usize,
        collapse_mode: PromoterWindowCollapseMode,
    ) -> Result<Vec<PromoterWindowRecord>, EngineError> {
        let records = {
            let dna = self
                .state
                .sequences
                .get(seq_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!("Sequence '{}' not found", seq_id),
                })?;
            self.derive_promoter_window_records(
                dna,
                gene_label,
                transcript_id,
                upstream_bp,
                downstream_bp,
                collapse_mode,
            )
        };
        let collapsed_records = Self::collapse_promoter_window_records_by_exact_span(records);
        if collapsed_records.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No transcript-derived promoter windows matched the requested filters on '{}'",
                    seq_id
                ),
            });
        }
        let dna = self
            .state
            .sequences
            .get_mut(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", seq_id),
            })?;
        let existing = dna
            .features()
            .iter()
            .filter_map(|feature| {
                let gene_label = Self::feature_qualifier_text(feature, "gene")
                    .or_else(|| Self::feature_qualifier_text(feature, "label"))
                    .unwrap_or_default();
                let generated = Self::feature_qualifier_text(feature, "gentle_generated")
                    .unwrap_or_default()
                    .trim()
                    .to_ascii_lowercase();
                let mut ranges = vec![];
                collect_location_ranges_usize(&feature.location, &mut ranges);
                let (start, end) = (
                    ranges.iter().map(|(start, _)| *start).min()?,
                    ranges.iter().map(|(_, end)| *end).max()?,
                );
                Some(format!(
                    "{}:{}:{}:{}:{}",
                    feature.kind.to_string().to_ascii_lowercase(),
                    generated,
                    if feature_is_reverse(feature) {
                        "-"
                    } else {
                        "+"
                    },
                    gene_label.to_ascii_lowercase(),
                    format!("{start}-{end}")
                ))
            })
            .collect::<HashSet<_>>();
        let mut added = 0usize;
        for record in &collapsed_records {
            let signature = format!(
                "promoter:{}:{}:{}:{}",
                ANNOTATE_PROMOTER_WINDOWS_GENERATED_TAG,
                record.strand,
                record
                    .gene_label
                    .clone()
                    .unwrap_or_else(|| record.transcript_label.clone())
                    .to_ascii_lowercase(),
                format!("{}-{}", record.start_0based, record.end_0based_exclusive)
            );
            if existing.contains(&signature) {
                continue;
            }
            dna.features_mut()
                .push(Self::build_promoter_window_feature(record));
            added += 1;
        }
        if added > 0 {
            Self::prepare_sequence(dna);
        }
        Ok(collapsed_records)
    }

    pub(crate) fn build_construct_reasoning_generated_promoter_evidence(
        &self,
        seq_id: &str,
        dna: &DNAsequence,
        upstream_bp: usize,
        downstream_bp: usize,
    ) -> Vec<DesignEvidence> {
        if dna.features().iter().any(|feature| {
            Self::construct_reasoning_role_from_feature(feature) == Some(ConstructRole::Promoter)
        }) {
            return vec![];
        }
        Self::collapse_promoter_window_records_by_exact_span(self.derive_promoter_window_records(
            dna,
            None,
            None,
            upstream_bp,
            downstream_bp,
            PromoterWindowCollapseMode::Transcript,
        ))
        .into_iter()
        .map(|record| DesignEvidence {
            evidence_id: format!(
                "generated_promoter_{}_{}_{}",
                Self::normalize_id_token(
                    record
                        .gene_label
                        .as_deref()
                        .unwrap_or(record.transcript_id.as_str())
                ),
                record.start_0based,
                record.end_0based_exclusive
            ),
            seq_id: seq_id.to_string(),
            start_0based: record.start_0based,
            end_0based_exclusive: record.end_0based_exclusive,
            strand: Some(record.strand.clone()),
            role: ConstructRole::Promoter,
            evidence_class: EvidenceClass::ContextEvidence,
            label: {
                let base_label = record
                    .gene_label
                    .clone()
                    .unwrap_or_else(|| record.transcript_label.clone());
                let transcript_count =
                    record.transcript_count.max(record.transcript_ids.len().max(1));
                if transcript_count > 1 {
                    format!("{base_label} promoter window ({transcript_count} tx)")
                } else {
                    base_label
                }
            },
            rationale: {
                let transcript_count =
                    record.transcript_count.max(record.transcript_ids.len().max(1));
                if transcript_count > 1 {
                    format!(
                        "Promoter window derived from transcript TSS geometry ({} bp upstream / {} bp downstream) and collapsed from {transcript_count} transcript interpretations because the sequence lacks imported promoter annotations.",
                        record.upstream_bp, record.downstream_bp
                    )
                } else {
                    format!(
                        "Promoter window derived from transcript TSS geometry ({} bp upstream / {} bp downstream) because the sequence lacks imported promoter annotations.",
                        record.upstream_bp, record.downstream_bp
                    )
                }
            },
            confidence: Some(0.7),
            context_tags: vec![
                "promoter".to_string(),
                "generated".to_string(),
                ANNOTATE_PROMOTER_WINDOWS_GENERATED_TAG.to_string(),
                "transcript_tss".to_string(),
            ],
            provenance_kind: "derived_promoter_window".to_string(),
            provenance_refs: record.transcript_ids.clone(),
            notes: vec![
                format!("tss_local_0based={}", record.tss_local_0based),
                format!(
                    "transcript_count={}",
                    record.transcript_count.max(record.transcript_ids.len().max(1))
                ),
                format!("transcript_ids={}", record.transcript_ids.join(",")),
                format!("transcript_labels={}", record.transcript_labels.join(",")),
            ],
            ..DesignEvidence::default()
        })
        .collect()
    }

    fn select_variant_feature<'a>(
        dna: &'a DNAsequence,
        variant_label_or_id: Option<&str>,
    ) -> Result<(usize, &'a gb_io::seq::Feature), EngineError> {
        let mut candidates = dna
            .features()
            .iter()
            .enumerate()
            .filter(|(_, feature)| {
                Self::construct_reasoning_role_from_feature(feature) == Some(ConstructRole::Variant)
            })
            .collect::<Vec<_>>();
        if candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: "No variant/variation feature is available on this sequence".to_string(),
            });
        }
        if let Some(query) = variant_label_or_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            candidates.retain(|(feature_id, feature)| {
                Self::feature_matches_identifier(
                    feature,
                    *feature_id,
                    query,
                    &["db_xref", "rs_id", "variant_id"],
                )
            });
            if candidates.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "No variant/variation feature on this sequence matches '{}'",
                        query
                    ),
                });
            }
        }
        candidates.sort_by(|left, right| {
            let left_span = Self::feature_span_bounds(left.1).unwrap_or((usize::MAX, usize::MAX));
            let right_span = Self::feature_span_bounds(right.1).unwrap_or((usize::MAX, usize::MAX));
            left_span
                .0
                .cmp(&right_span.0)
                .then_with(|| left_span.1.cmp(&right_span.1))
                .then_with(|| {
                    Self::feature_display_label(left.1, left.0)
                        .to_ascii_lowercase()
                        .cmp(&Self::feature_display_label(right.1, right.0).to_ascii_lowercase())
                })
        });
        candidates.into_iter().next().ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: "No variant/variation feature is available on this sequence".to_string(),
        })
    }

    fn promoter_window_signed_tss_distance(
        record: &PromoterWindowRecord,
        variant_start_0based: usize,
    ) -> isize {
        if record.strand == "-" {
            record.tss_local_0based as isize - variant_start_0based as isize
        } else {
            variant_start_0based as isize - record.tss_local_0based as isize
        }
    }

    fn promoter_window_overlaps_variant(
        record: &PromoterWindowRecord,
        variant_start_0based: usize,
        variant_end_0based_exclusive: usize,
    ) -> bool {
        record.end_0based_exclusive > variant_start_0based
            && record.start_0based < variant_end_0based_exclusive
    }

    fn choose_promoter_context_record<'a>(
        records: &'a [PromoterWindowRecord],
        variant_start_0based: usize,
        variant_end_0based_exclusive: usize,
        transcript_filter_present: bool,
    ) -> (
        Option<&'a PromoterWindowRecord>,
        String,
        Option<isize>,
        bool,
    ) {
        if records.is_empty() {
            return (None, "no_transcript_context".to_string(), None, false);
        }
        let mut indexed = records
            .iter()
            .map(|record| {
                let overlap = Self::promoter_window_overlaps_variant(
                    record,
                    variant_start_0based,
                    variant_end_0based_exclusive,
                );
                let signed_distance =
                    Self::promoter_window_signed_tss_distance(record, variant_start_0based);
                (record, overlap, signed_distance)
            })
            .collect::<Vec<_>>();
        indexed.sort_by(|left, right| {
            right
                .1
                .cmp(&left.1)
                .then_with(|| left.2.abs().cmp(&right.2.abs()))
                .then_with(|| {
                    left.0
                        .gene_label
                        .clone()
                        .unwrap_or_else(|| left.0.transcript_label.clone())
                        .to_ascii_lowercase()
                        .cmp(
                            &right
                                .0
                                .gene_label
                                .clone()
                                .unwrap_or_else(|| right.0.transcript_label.clone())
                                .to_ascii_lowercase(),
                        )
                })
                .then_with(|| {
                    left.0
                        .transcript_id
                        .to_ascii_lowercase()
                        .cmp(&right.0.transcript_id.to_ascii_lowercase())
                })
        });
        let best = indexed.first().copied();
        let transcript_ambiguity_status = match (records.len(), transcript_filter_present, best) {
            (_, _, None) => "no_transcript_context".to_string(),
            (1, true, Some(_)) => "transcript_filtered".to_string(),
            (1, false, Some(_)) => "single_transcript".to_string(),
            (_, _, Some(best)) => {
                let tied_count = indexed
                    .iter()
                    .filter(|row| row.1 == best.1 && row.2.abs() == best.2.abs())
                    .count();
                if tied_count > 1 {
                    "multi_transcript_ambiguous".to_string()
                } else {
                    "multi_transcript_defaulted".to_string()
                }
            }
        };
        best.map_or(
            (None, transcript_ambiguity_status.clone(), None, false),
            |row| (Some(row.0), transcript_ambiguity_status, Some(row.2), row.1),
        )
    }

    fn variant_promoter_overlap_rows(
        dna: &DNAsequence,
        promoter_windows: &[PromoterWindowRecord],
        variant_start_0based: usize,
        variant_end_0based_exclusive: usize,
    ) -> Vec<VariantPromoterContextEvidenceRow> {
        let mut rows = vec![];
        let mut seen = HashSet::new();
        for (feature_id, feature) in dna.features().iter().enumerate() {
            let Some(role) = Self::construct_reasoning_role_from_feature(feature) else {
                continue;
            };
            if !matches!(
                role,
                ConstructRole::Gene
                    | ConstructRole::Transcript
                    | ConstructRole::Promoter
                    | ConstructRole::Tfbs
            ) || !Self::feature_overlaps_span(
                feature,
                variant_start_0based,
                variant_end_0based_exclusive,
            ) {
                continue;
            }
            let Some((start_0based, end_0based_exclusive)) = Self::feature_span_bounds(feature)
            else {
                continue;
            };
            let source = Self::feature_qualifier_text(feature, "gentle_generated")
                .filter(|value| !value.trim().is_empty())
                .unwrap_or_else(|| "annotation".to_string());
            let row = VariantPromoterContextEvidenceRow {
                role: role.as_str().to_string(),
                kind: feature.kind.to_string(),
                label: Self::feature_display_label(feature, feature_id),
                start_0based,
                end_0based_exclusive,
                strand: Some(
                    if feature_is_reverse(feature) {
                        "-"
                    } else {
                        "+"
                    }
                    .to_string(),
                ),
                source,
            };
            let key = format!(
                "{}:{}:{}:{}:{}",
                row.role, row.kind, row.label, row.start_0based, row.end_0based_exclusive
            );
            if seen.insert(key) {
                rows.push(row);
            }
        }
        for record in promoter_windows.iter().filter(|record| {
            Self::promoter_window_overlaps_variant(
                record,
                variant_start_0based,
                variant_end_0based_exclusive,
            )
        }) {
            let label = record
                .gene_label
                .clone()
                .unwrap_or_else(|| record.transcript_label.clone());
            let row = VariantPromoterContextEvidenceRow {
                role: ConstructRole::Promoter.as_str().to_string(),
                kind: "promoter".to_string(),
                label,
                start_0based: record.start_0based,
                end_0based_exclusive: record.end_0based_exclusive,
                strand: Some(record.strand.clone()),
                source: ANNOTATE_PROMOTER_WINDOWS_GENERATED_TAG.to_string(),
            };
            let key = format!(
                "{}:{}:{}:{}:{}",
                row.role, row.kind, row.label, row.start_0based, row.end_0based_exclusive
            );
            if seen.insert(key) {
                rows.push(row);
            }
        }
        rows.sort_by(|left, right| {
            left.start_0based
                .cmp(&right.start_0based)
                .then_with(|| left.end_0based_exclusive.cmp(&right.end_0based_exclusive))
                .then_with(|| left.role.cmp(&right.role))
                .then_with(|| {
                    left.label
                        .to_ascii_lowercase()
                        .cmp(&right.label.to_ascii_lowercase())
                })
        });
        rows
    }

    pub(crate) fn summarize_variant_promoter_context(
        &self,
        input: &str,
        variant_label_or_id: Option<&str>,
        gene_label: Option<&str>,
        transcript_id: Option<&str>,
        promoter_upstream_bp: usize,
        promoter_downstream_bp: usize,
        tfbs_focus_half_window_bp: usize,
    ) -> Result<VariantPromoterContextReport, EngineError> {
        let dna = self.state.sequences.get(input).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!("Sequence '{}' not found", input),
        })?;
        let (variant_feature_id, variant_feature) =
            Self::select_variant_feature(dna, variant_label_or_id)?;
        let (variant_start_0based, variant_end_0based_exclusive) =
            Self::feature_span_bounds(variant_feature).ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: "Selected variant feature has no usable coordinates".to_string(),
            })?;
        let variant_label = Self::feature_display_label(variant_feature, variant_feature_id);
        let promoter_windows = self.derive_promoter_window_records(
            dna,
            gene_label,
            transcript_id,
            promoter_upstream_bp,
            promoter_downstream_bp,
            PromoterWindowCollapseMode::Transcript,
        );
        let (chosen_record, transcript_ambiguity_status, signed_tss_distance_bp, promoter_overlap) =
            Self::choose_promoter_context_record(
                &promoter_windows,
                variant_start_0based,
                variant_end_0based_exclusive,
                transcript_id.is_some(),
            );

        let mut evidence = self.build_construct_reasoning_sequence_evidence_with_promoter_params(
            input,
            dna,
            promoter_upstream_bp,
            promoter_downstream_bp,
        );
        if !evidence.iter().any(|row| {
            row.role == ConstructRole::Promoter
                && row.evidence_id.starts_with("generated_promoter_")
        }) && !dna.features().iter().any(|feature| {
            Self::construct_reasoning_role_from_feature(feature) == Some(ConstructRole::Promoter)
        }) {
            evidence.extend(self.build_construct_reasoning_generated_promoter_evidence(
                input,
                dna,
                promoter_upstream_bp,
                promoter_downstream_bp,
            ));
        }
        let variant_evidence = evidence
            .iter()
            .find(|row| {
                row.role == ConstructRole::Variant
                    && row.start_0based == variant_start_0based
                    && row.end_0based_exclusive == variant_end_0based_exclusive
            })
            .cloned()
            .unwrap_or_else(|| DesignEvidence {
                evidence_id: format!(
                    "variant_{}_{}_{}",
                    Self::normalize_id_token(&variant_label),
                    variant_start_0based,
                    variant_end_0based_exclusive
                ),
                seq_id: input.to_string(),
                start_0based: variant_start_0based,
                end_0based_exclusive: variant_end_0based_exclusive,
                strand: Some(
                    if feature_is_reverse(variant_feature) {
                        "-"
                    } else {
                        "+"
                    }
                    .to_string(),
                ),
                role: ConstructRole::Variant,
                evidence_class: EvidenceClass::ReliableAnnotation,
                label: variant_label.clone(),
                rationale: "Variant marker selected for promoter-context summarization."
                    .to_string(),
                ..DesignEvidence::default()
            });
        let variant_summary = Self::construct_reasoning_variant_summary(
            dna,
            variant_feature,
            &variant_evidence,
            &evidence,
        );

        let tfbs_present = dna.features().iter().any(Self::is_tfbs_feature);
        let tfbs_region_summary = if tfbs_present {
            let focus_start_0based = variant_start_0based.saturating_sub(tfbs_focus_half_window_bp);
            let focus_end_0based_exclusive = variant_end_0based_exclusive
                .saturating_add(tfbs_focus_half_window_bp)
                .min(dna.len());
            Some(self.summarize_tfbs_region(TfbsRegionSummaryRequest {
                seq_id: input.to_string(),
                focus_start_0based,
                focus_end_0based_exclusive:
                    focus_end_0based_exclusive.max(focus_start_0based.saturating_add(1)),
                context_start_0based: Some(0),
                context_end_0based_exclusive: Some(dna.len()),
                min_focus_occurrences: 1,
                min_context_occurrences: 0,
                limit: Some(TFBS_REGION_SUMMARY_DEFAULT_LIMIT.min(25)),
            })?)
        } else {
            None
        };
        let tfbs_near_variant_status = match tfbs_region_summary.as_ref() {
            Some(summary) if summary.focus_hit_count > 0 => "tfbs_annotations_near_variant",
            Some(_) => "tfbs_annotations_absent_in_focus_window",
            None => "no_tfbs_annotations_present",
        }
        .to_string();

        let overlapping_evidence = Self::variant_promoter_overlap_rows(
            dna,
            &promoter_windows,
            variant_start_0based,
            variant_end_0based_exclusive,
        );
        let overlapping_gene_labels = overlapping_evidence
            .iter()
            .filter(|row| row.role == ConstructRole::Gene.as_str())
            .map(|row| row.label.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let overlapping_transcript_labels = overlapping_evidence
            .iter()
            .filter(|row| row.role == ConstructRole::Transcript.as_str())
            .map(|row| row.label.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let overlapping_promoter_labels = overlapping_evidence
            .iter()
            .filter(|row| row.role == ConstructRole::Promoter.as_str())
            .map(|row| row.label.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let overlapping_tfbs_labels = tfbs_region_summary
            .as_ref()
            .map(|summary| {
                summary
                    .rows
                    .iter()
                    .map(|row| row.tf_name.clone())
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        let rationale = if let Some(record) = chosen_record {
            let promoter_text = if promoter_overlap {
                "falls inside"
            } else {
                "does not overlap"
            };
            format!(
                "Variant '{}' {} the transcript-derived promoter window for '{}' (transcript '{}', signed_tss_distance_bp={}). GENtle uses this as a deterministic promoter-context handoff for reporter-design planning rather than as wet-lab validation.",
                variant_label,
                promoter_text,
                record
                    .gene_label
                    .clone()
                    .unwrap_or_else(|| record.transcript_label.clone()),
                record.transcript_id,
                signed_tss_distance_bp.unwrap_or_default()
            )
        } else {
            format!(
                "Variant '{}' was located on '{}', but no transcript-derived promoter window matched the requested filters. GENtle therefore records the locus context without claiming promoter overlap.",
                variant_label, input
            )
        };

        Ok(VariantPromoterContextReport {
            schema: VARIANT_PROMOTER_CONTEXT_SCHEMA.to_string(),
            seq_id: input.to_string(),
            sequence_length_bp: dna.len(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            variant_label,
            variant_feature_id: Some(variant_feature_id),
            variant_start_0based,
            variant_end_0based_exclusive,
            variant_class: Self::feature_qualifier_text(variant_feature, "vcf_variant_class"),
            genomic_ref: Self::feature_qualifier_text(variant_feature, "vcf_ref"),
            genomic_alt: Self::feature_qualifier_text(variant_feature, "vcf_alt"),
            genome_anchor: self.sequence_genome_anchor_summary(input).ok(),
            chosen_gene_label: chosen_record.and_then(|record| record.gene_label.clone()),
            chosen_transcript_id: chosen_record.map(|record| record.transcript_id.clone()),
            chosen_transcript_label: chosen_record.map(|record| record.transcript_label.clone()),
            transcript_ambiguity_status,
            promoter_upstream_bp,
            promoter_downstream_bp,
            promoter_overlap,
            signed_tss_distance_bp,
            overlapping_gene_labels,
            overlapping_transcript_labels,
            overlapping_promoter_labels,
            overlapping_tfbs_labels,
            overlapping_evidence,
            promoter_windows_considered: promoter_windows,
            effect_tags: variant_summary.effect_tags,
            suggested_assay_ids: variant_summary.suggested_assay_ids,
            tfbs_focus_half_window_bp,
            tfbs_near_variant_status,
            tfbs_region_summary,
            rationale,
        })
    }

    pub(crate) fn write_variant_promoter_context_json(
        &self,
        report: &VariantPromoterContextReport,
        path: &str,
    ) -> Result<(), EngineError> {
        let text = serde_json::to_string_pretty(report).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize variant-promoter context report '{}' for '{}': {e}",
                report.seq_id, path
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write variant-promoter context report to '{path}': {e}"),
        })
    }

    pub(crate) fn suggest_promoter_reporter_fragments(
        &self,
        input: &str,
        variant_label_or_id: Option<&str>,
        gene_label: Option<&str>,
        transcript_id: Option<&str>,
        retain_downstream_from_tss_bp: usize,
        retain_upstream_beyond_variant_bp: usize,
        max_candidates: usize,
    ) -> Result<PromoterReporterCandidateSet, EngineError> {
        if max_candidates == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "max_candidates must be >= 1".to_string(),
            });
        }
        let context = self.summarize_variant_promoter_context(
            input,
            variant_label_or_id,
            gene_label,
            transcript_id,
            DEFAULT_PROMOTER_WINDOW_UPSTREAM_BP,
            DEFAULT_PROMOTER_WINDOW_DOWNSTREAM_BP,
            DEFAULT_VARIANT_PROMOTER_TFBS_FOCUS_HALF_WINDOW_BP,
        )?;
        let mut candidates = context
            .promoter_windows_considered
            .iter()
            .map(|record| {
                let signed_tss_distance_bp =
                    Self::promoter_window_signed_tss_distance(record, context.variant_start_0based);
                let promoter_overlap = Self::promoter_window_overlaps_variant(
                    record,
                    context.variant_start_0based,
                    context.variant_end_0based_exclusive,
                );
                let (start_0based, end_0based_exclusive) = if record.strand == "-" {
                    (
                        record
                            .tss_local_0based
                            .saturating_sub(retain_downstream_from_tss_bp),
                        context
                            .variant_end_0based_exclusive
                            .saturating_add(retain_upstream_beyond_variant_bp)
                            .min(context.sequence_length_bp),
                    )
                } else {
                    (
                        context
                            .variant_start_0based
                            .saturating_sub(retain_upstream_beyond_variant_bp),
                        record
                            .tss_local_0based
                            .saturating_add(retain_downstream_from_tss_bp)
                            .saturating_add(1)
                            .min(context.sequence_length_bp),
                    )
                };
                let length_bp = end_0based_exclusive.saturating_sub(start_0based);
                let candidate_id = format!(
                    "{}_{}_{}_{}_promoter_fragment",
                    Self::normalize_id_token(input),
                    Self::normalize_id_token(&record.transcript_id),
                    start_0based.saturating_add(1),
                    end_0based_exclusive
                );
                let rationale = if promoter_overlap {
                    format!(
                        "Keeps {} bp downstream of the TSS and {} bp beyond the SNP on the biologically upstream side while retaining promoter overlap for transcript '{}'.",
                        retain_downstream_from_tss_bp,
                        retain_upstream_beyond_variant_bp,
                        record.transcript_id
                    )
                } else {
                    format!(
                        "Transcript '{}' is the nearest deterministic promoter-context candidate even though the SNP does not lie inside the derived promoter window.",
                        record.transcript_id
                    )
                };
                PromoterReporterFragmentCandidate {
                    candidate_id,
                    gene_label: record.gene_label.clone(),
                    transcript_id: record.transcript_id.clone(),
                    transcript_label: record.transcript_label.clone(),
                    strand: record.strand.clone(),
                    tss_local_0based: record.tss_local_0based,
                    variant_start_0based: context.variant_start_0based,
                    variant_end_0based_exclusive: context.variant_end_0based_exclusive,
                    start_0based,
                    end_0based_exclusive,
                    length_bp,
                    retain_downstream_from_tss_bp,
                    retain_upstream_beyond_variant_bp,
                    promoter_overlap,
                    signed_tss_distance_bp,
                    rank: 0,
                    recommended: false,
                    rationale,
                }
            })
            .filter(|row| row.end_0based_exclusive > row.start_0based)
            .collect::<Vec<_>>();
        if candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No promoter-reporter fragment candidates could be derived for '{}'",
                    input
                ),
            });
        }
        candidates.sort_by(|left, right| {
            right
                .promoter_overlap
                .cmp(&left.promoter_overlap)
                .then_with(|| {
                    left.signed_tss_distance_bp
                        .abs()
                        .cmp(&right.signed_tss_distance_bp.abs())
                })
                .then_with(|| {
                    left.gene_label
                        .clone()
                        .unwrap_or_else(|| left.transcript_label.clone())
                        .to_ascii_lowercase()
                        .cmp(
                            &right
                                .gene_label
                                .clone()
                                .unwrap_or_else(|| right.transcript_label.clone())
                                .to_ascii_lowercase(),
                        )
                })
                .then_with(|| {
                    left.transcript_id
                        .to_ascii_lowercase()
                        .cmp(&right.transcript_id.to_ascii_lowercase())
                })
        });
        candidates.truncate(max_candidates);
        for (idx, candidate) in candidates.iter_mut().enumerate() {
            candidate.rank = idx + 1;
            candidate.recommended = idx == 0;
        }
        let recommended = candidates.first().cloned().ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "No promoter-reporter fragment candidates could be derived for '{}'",
                input
            ),
        })?;
        Ok(PromoterReporterCandidateSet {
            schema: PROMOTER_REPORTER_CANDIDATES_SCHEMA.to_string(),
            seq_id: input.to_string(),
            sequence_length_bp: context.sequence_length_bp,
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            variant_label: context.variant_label,
            chosen_gene_label: recommended.gene_label.clone(),
            chosen_transcript_id: Some(recommended.transcript_id.clone()),
            chosen_transcript_label: Some(recommended.transcript_label.clone()),
            transcript_ambiguity_status: context.transcript_ambiguity_status,
            retain_downstream_from_tss_bp,
            retain_upstream_beyond_variant_bp,
            max_candidates,
            recommended_candidate_id: recommended.candidate_id.clone(),
            suggested_assay_ids: context.suggested_assay_ids,
            candidates,
        })
    }

    pub(crate) fn write_promoter_reporter_candidates_json(
        &self,
        report: &PromoterReporterCandidateSet,
        path: &str,
    ) -> Result<(), EngineError> {
        let text = serde_json::to_string_pretty(report).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize promoter-reporter candidates '{}' for '{}': {e}",
                report.seq_id, path
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write promoter-reporter candidate report to '{path}': {e}"),
        })
    }

    fn parse_single_base_alternate_allele(raw_alt: &str) -> Result<String, EngineError> {
        let trimmed = raw_alt.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Selected variant feature has no alternate allele (vcf_alt)".to_string(),
            });
        }
        let alts = trimmed
            .split([',', '/', '|'])
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .collect::<Vec<_>>();
        if alts.len() != 1 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "MaterializeVariantAllele currently supports only single alternate alleles; observed '{}'",
                    raw_alt
                ),
            });
        }
        let alt = alts[0].to_ascii_uppercase();
        if alt.len() != 1 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "MaterializeVariantAllele currently supports only SNVs; alternate allele '{}' is not one base",
                    alt
                ),
            });
        }
        Ok(alt)
    }

    pub(crate) fn materialize_variant_allele_sequence(
        &self,
        input: &str,
        variant_label_or_id: Option<&str>,
        allele: VariantAlleleChoice,
        output_id: Option<&str>,
    ) -> Result<(String, DNAsequence), EngineError> {
        let dna = self.state.sequences.get(input).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!("Sequence '{}' not found", input),
        })?;
        let (variant_feature_id, variant_feature) =
            Self::select_variant_feature(dna, variant_label_or_id)?;
        let (variant_start_0based, variant_end_0based_exclusive) =
            Self::feature_span_bounds(variant_feature).ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: "Selected variant feature has no usable coordinates".to_string(),
            })?;
        if variant_end_0based_exclusive.saturating_sub(variant_start_0based) != 1 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "MaterializeVariantAllele currently supports only single-nucleotide variants"
                        .to_string(),
            });
        }
        let reference = Self::feature_qualifier_text(variant_feature, "vcf_ref")
            .map(|value| value.trim().to_ascii_uppercase())
            .filter(|value| value.len() == 1)
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "Selected variant feature does not carry a single-base reference allele (vcf_ref)"
                        .to_string(),
            })?;
        let alternate_raw =
            Self::feature_qualifier_text(variant_feature, "vcf_alt").ok_or_else(|| {
                EngineError {
                    code: ErrorCode::InvalidInput,
                    message:
                        "Selected variant feature does not carry an alternate allele (vcf_alt)"
                            .to_string(),
                }
            })?;
        let alternate = Self::parse_single_base_alternate_allele(&alternate_raw)?;
        let current_base = dna
            .forward_bytes()
            .get(variant_start_0based)
            .copied()
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: "Selected variant position falls outside the sequence".to_string(),
            })? as char;
        let reference_base = reference.chars().next().unwrap_or('N');
        if !current_base.eq_ignore_ascii_case(&reference_base) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Sequence base '{}' at position {} does not match vcf_ref '{}'",
                    current_base,
                    variant_start_0based.saturating_add(1),
                    reference_base
                ),
            });
        }
        let materialized_base = match allele {
            VariantAlleleChoice::Reference => reference.clone(),
            VariantAlleleChoice::Alternate => alternate.clone(),
        };
        let mut seq = dna.clone_seq_record();
        if let Some(base) = materialized_base.as_bytes().first() {
            seq.seq[variant_start_0based] = *base;
        }
        if let Some(feature) = seq.features.get_mut(variant_feature_id) {
            feature.qualifiers.retain(|(key, _)| {
                let key_text = key.to_string();
                !matches!(
                    key_text.as_str(),
                    "materialized_allele" | "materialized_base" | "materialized_from_variant"
                )
            });
            feature.qualifiers.push((
                "materialized_allele".into(),
                Some(allele.as_str().to_string()),
            ));
            feature
                .qualifiers
                .push(("materialized_base".into(), Some(materialized_base.clone())));
            feature.qualifiers.push((
                "materialized_from_variant".into(),
                Some(Self::feature_display_label(
                    variant_feature,
                    variant_feature_id,
                )),
            ));
        }
        let mut out = DNAsequence::from_genbank_seq(seq);
        Self::prepare_sequence(&mut out);
        let variant_token = Self::normalize_id_token(&Self::feature_display_label(
            variant_feature,
            variant_feature_id,
        ));
        let base = output_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string())
            .unwrap_or_else(|| {
                format!(
                    "{}_{}_{}",
                    Self::local_strip_dbsnp_auto_region_id_to_rsid(input)
                        .unwrap_or_else(|| input.to_string()),
                    variant_token,
                    allele.as_str()
                )
            });
        Ok((base, out))
    }

    fn local_strip_dbsnp_auto_region_id_to_rsid(input: &str) -> Option<String> {
        let mut parts = input.split('_').collect::<Vec<_>>();
        if parts.len() < 4 {
            return None;
        }
        let rs_id = parts.first().copied()?.trim();
        if !(rs_id.starts_with("rs") && rs_id[2..].chars().all(|ch| ch.is_ascii_digit())) {
            return None;
        }
        let last = parts.pop()?.trim();
        let penultimate = parts.pop()?.trim();
        if !last.chars().all(|ch| ch.is_ascii_digit())
            || !penultimate.chars().all(|ch| ch.is_ascii_digit())
        {
            return None;
        }
        Some(rs_id.to_string())
    }
}
