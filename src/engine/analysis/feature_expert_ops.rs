//! Feature-expert, splicing, and isoform helper routines used by engine
//! operations.
//!
//! This module keeps feature-oriented interpretation logic close together so
//! inspection/render/export operations reuse the same qualifier parsing rules.

use super::*;

impl GentleEngine {
    pub(super) fn feature_qualifier_text(
        feature: &gb_io::seq::Feature,
        key: &str,
    ) -> Option<String> {
        feature
            .qualifier_values(key.into())
            .map(|value| value.split_whitespace().collect::<Vec<_>>().join(" "))
            .map(|value| value.trim().to_string())
            .find(|value| !value.is_empty())
    }

    pub(super) fn feature_qualifier_f64(feature: &gb_io::seq::Feature, key: &str) -> Option<f64> {
        feature
            .qualifier_values(key.into())
            .next()
            .and_then(|v| v.trim().parse::<f64>().ok())
    }

    pub(super) fn first_nonempty_feature_qualifier(
        feature: &gb_io::seq::Feature,
        keys: &[&str],
    ) -> Option<String> {
        for key in keys {
            if let Some(value) = Self::feature_qualifier_text(feature, key) {
                return Some(value);
            }
        }
        None
    }

    pub(super) fn serialize_ranges_1based(ranges: &[(usize, usize)]) -> Option<String> {
        let parts = ranges
            .iter()
            .filter_map(|(start, end)| {
                (*start > 0 && *end >= *start).then_some(format!("{start}-{end}"))
            })
            .collect::<Vec<_>>();
        (!parts.is_empty()).then_some(parts.join(","))
    }

    pub(super) fn parse_ranges_1based(raw: &str) -> Vec<(usize, usize)> {
        let mut out = vec![];
        for token in raw.split(',') {
            let trimmed = token.trim();
            if trimmed.is_empty() {
                continue;
            }
            let mut pieces = trimmed.splitn(2, '-');
            let Some(start_raw) = pieces.next() else {
                continue;
            };
            let Some(end_raw) = pieces.next() else {
                continue;
            };
            let Ok(start) = start_raw.trim().parse::<usize>() else {
                continue;
            };
            let Ok(end) = end_raw.trim().parse::<usize>() else {
                continue;
            };
            if start == 0 || end < start {
                continue;
            }
            out.push((start, end));
        }
        out.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        out.dedup();
        out
    }

    pub(super) fn feature_qualifier_ranges_0based(
        feature: &gb_io::seq::Feature,
        key: &str,
    ) -> Vec<(usize, usize)> {
        Self::feature_qualifier_text(feature, key)
            .map(|raw| {
                Self::parse_ranges_1based(&raw)
                    .into_iter()
                    .map(|(start_1based, end_1based)| (start_1based - 1, end_1based))
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default()
    }

    pub(super) fn is_tfbs_feature(feature: &gb_io::seq::Feature) -> bool {
        matches!(
            feature.kind.to_string().to_ascii_uppercase().as_str(),
            "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND"
        )
    }

    pub(super) fn feature_display_label(
        feature: &gb_io::seq::Feature,
        feature_id: usize,
    ) -> String {
        let fallback = format!("{} #{}", feature.kind.to_string(), feature_id + 1);
        let label = Self::first_nonempty_feature_qualifier(
            feature,
            &[
                "label",
                "name",
                "standard_name",
                "gene",
                "protein_id",
                "product",
                "region_name",
                "bound_moiety",
            ],
        )
        .unwrap_or(fallback);
        let trimmed = label.trim();
        if trimmed.is_empty() {
            format!("{} #{}", feature.kind.to_string(), feature_id + 1)
        } else {
            trimmed.to_string()
        }
    }

    pub(super) fn normalize_column_frequencies(counts: [f64; 4]) -> [f64; 4] {
        let total = counts
            .iter()
            .copied()
            .map(|v| if v.is_finite() && v > 0.0 { v } else { 0.0 })
            .sum::<f64>();
        if total <= 0.0 {
            [0.25_f64, 0.25_f64, 0.25_f64, 0.25_f64]
        } else {
            [
                (counts[0].max(0.0) / total).clamp(0.0, 1.0),
                (counts[1].max(0.0) / total).clamp(0.0, 1.0),
                (counts[2].max(0.0) / total).clamp(0.0, 1.0),
                (counts[3].max(0.0) / total).clamp(0.0, 1.0),
            ]
        }
    }

    pub(super) fn information_content_bits(frequencies: &[f64; 4]) -> f64 {
        let entropy = frequencies
            .iter()
            .copied()
            .filter(|p| *p > 0.0)
            .map(|p| -p * p.log2())
            .sum::<f64>();
        (2.0 - entropy).clamp(0.0, 2.0)
    }

    pub(super) fn resolve_tfbs_scoring_motif(
        feature: &gb_io::seq::Feature,
    ) -> Result<(String, Option<String>, Vec<[f64; 4]>), EngineError> {
        let label = Self::feature_qualifier_text(feature, "label");
        let mut tokens = vec![];
        if let Some(tf_id) = Self::feature_qualifier_text(feature, "tf_id") {
            tokens.push(tf_id);
        }
        if let Some(bound) = Self::first_nonempty_feature_qualifier(
            feature,
            &["bound_moiety", "standard_name", "gene", "name"],
        ) {
            tokens.push(bound);
        }
        if let Some(raw_label) = label {
            let trimmed = raw_label.trim().to_string();
            if !trimmed.is_empty() {
                tokens.push(trimmed.clone());
                let label_upper = trimmed.to_ascii_uppercase();
                if let Some(stripped) = label_upper.strip_prefix("TFBS ") {
                    tokens.push(stripped.trim().to_string());
                }
            }
        }

        let mut last_error: Option<EngineError> = None;
        for token in tokens {
            if token.trim().is_empty() {
                continue;
            }
            match Self::resolve_tf_motif_for_scoring(&token) {
                Ok((tf_id, tf_name, _, matrix_counts)) => {
                    return Ok((tf_id, tf_name, matrix_counts));
                }
                Err(err) => {
                    last_error = Some(err);
                }
            }
        }
        Err(last_error.unwrap_or(EngineError {
            code: ErrorCode::InvalidInput,
            message:
                "Could not resolve TF motif for feature (missing tf_id/label resolvable in motif registry)"
                    .to_string(),
        }))
    }

    pub(super) fn build_tfbs_expert_view(
        &self,
        seq_id: &str,
        feature_id: usize,
    ) -> Result<TfbsExpertView, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let feature = dna.features().get(feature_id).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "Feature id '{}' was not found in sequence '{}'",
                feature_id, seq_id
            ),
        })?;
        if !Self::is_tfbs_feature(feature) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Feature '{}' in '{}' is not a TFBS/protein-bind feature",
                    feature_id, seq_id
                ),
            });
        }
        let (from, _to) = feature.location.find_bounds().map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not read TFBS feature range: {e}"),
        })?;
        if from < 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "TFBS feature has negative range bounds".to_string(),
            });
        }

        let (tf_id, tf_name, matrix_counts) = Self::resolve_tfbs_scoring_motif(feature)?;
        if matrix_counts.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Resolved motif '{tf_id}' has empty matrix"),
            });
        }
        let motif_length = matrix_counts.len();
        let start = from as usize;
        let end = start.saturating_add(motif_length);
        let mut matched_bytes = dna.get_range_safe(start..end).ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "TFBS window {}..{} for feature '{}' exceeds sequence bounds",
                start + 1,
                end,
                feature_id
            ),
        })?;
        let is_reverse = feature_is_reverse(feature);
        if is_reverse {
            matched_bytes = Self::reverse_complement_bytes(&matched_bytes);
        }
        let matched_sequence = String::from_utf8_lossy(&matched_bytes).to_ascii_uppercase();

        let (llr_matrix, true_log_odds_matrix) = Self::prepare_scoring_matrices(&matrix_counts);
        let mut columns: Vec<TfbsExpertColumn> = Vec::with_capacity(motif_length);
        for idx in 0..motif_length {
            let counts = matrix_counts[idx];
            let frequencies = Self::normalize_column_frequencies(counts);
            let information_content_bits = Self::information_content_bits(&frequencies);
            let match_base = matched_bytes
                .get(idx)
                .map(|b| (*b as char).to_ascii_uppercase());
            let match_idx = matched_bytes.get(idx).and_then(|b| Self::base_to_idx(*b));
            let (match_frequency, llr_bits, true_log_odds_bits) = if let Some(base_idx) = match_idx
            {
                (
                    Some(frequencies[base_idx]),
                    Some(llr_matrix[idx][base_idx]),
                    Some(true_log_odds_matrix[idx][base_idx]),
                )
            } else {
                (None, None, None)
            };
            columns.push(TfbsExpertColumn {
                index_1based: idx + 1,
                counts,
                frequencies,
                information_content_bits,
                match_base,
                match_frequency,
                llr_bits,
                true_log_odds_bits,
                llr_rank_desc: None,
            });
        }

        let mut ranked = columns
            .iter()
            .enumerate()
            .filter_map(|(idx, col)| col.llr_bits.map(|score| (idx, score)))
            .collect::<Vec<_>>();
        ranked.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        for (rank_idx, (col_idx, _)) in ranked.iter().enumerate() {
            if let Some(col) = columns.get_mut(*col_idx) {
                col.llr_rank_desc = Some(rank_idx + 1);
            }
        }

        let llr_from_columns = columns.iter().filter_map(|col| col.llr_bits).sum::<f64>();
        let true_log_odds_from_columns = columns
            .iter()
            .filter_map(|col| col.true_log_odds_bits)
            .sum::<f64>();
        let llr_total_bits =
            Self::feature_qualifier_f64(feature, "llr_bits").or(Some(llr_from_columns));
        let llr_quantile = Self::feature_qualifier_f64(feature, "llr_quantile");
        let true_log_odds_total_bits = Self::feature_qualifier_f64(feature, "true_log_odds_bits")
            .or_else(|| Self::feature_qualifier_f64(feature, "log_odds_ratio_bits"))
            .or(Some(true_log_odds_from_columns));
        let true_log_odds_quantile = Self::feature_qualifier_f64(feature, "true_log_odds_quantile")
            .or_else(|| Self::feature_qualifier_f64(feature, "log_odds_ratio_quantile"));

        Ok(TfbsExpertView {
            seq_id: seq_id.to_string(),
            feature_id,
            feature_label: Self::feature_display_label(feature, feature_id),
            tf_id,
            tf_name,
            strand: if is_reverse {
                "-".to_string()
            } else {
                "+".to_string()
            },
            start_1based: start + 1,
            end_1based: end,
            motif_length,
            matched_sequence,
            llr_total_bits,
            llr_quantile,
            true_log_odds_total_bits,
            true_log_odds_quantile,
            instruction: TFBS_EXPERT_INSTRUCTION.to_string(),
            columns,
        })
    }

    pub(super) fn resolve_restriction_site_target(
        &self,
        seq_id: &str,
        cut_pos_1based: usize,
        enzyme: Option<&str>,
        recognition_start_1based: Option<usize>,
        recognition_end_1based: Option<usize>,
    ) -> Result<(RestrictionEnzymeKey, Vec<String>, Option<String>), EngineError> {
        if cut_pos_1based == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Restriction-site cut_pos_1based must be >= 1".to_string(),
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
        let enzyme_filter = enzyme
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_ascii_uppercase());
        let mut candidates: Vec<(RestrictionEnzymeKey, Vec<String>)> = Vec::new();
        for (key, names) in dna.restriction_enzyme_groups() {
            if key.pos() < 0 || (key.pos() as usize + 1) != cut_pos_1based {
                continue;
            }
            if let Some(start) = recognition_start_1based {
                if key.from() < 0 || (key.from() as usize + 1) != start {
                    continue;
                }
            }
            if let Some(end) = recognition_end_1based {
                if key.to() < 0 || key.to() as usize != end {
                    continue;
                }
            }
            if let Some(filter) = &enzyme_filter {
                if !names
                    .iter()
                    .any(|name| name.to_ascii_uppercase() == *filter)
                {
                    continue;
                }
            }
            candidates.push((key.clone(), names.clone()));
        }
        if candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No restriction-site match found for cut position {} in '{}'",
                    cut_pos_1based, seq_id
                ),
            });
        }
        if candidates.len() > 1 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Restriction-site target is ambiguous for cut position {} in '{}'; provide enzyme and/or recognition start/end",
                    cut_pos_1based, seq_id
                ),
            });
        }
        let (key, mut names) = candidates.into_iter().next().unwrap_or_default();
        names.sort_by(|a, b| a.to_ascii_uppercase().cmp(&b.to_ascii_uppercase()));
        names.dedup_by(|left, right| left.eq_ignore_ascii_case(right));
        let selected_enzyme = if let Some(filter) = enzyme_filter {
            names
                .iter()
                .find(|name| name.to_ascii_uppercase() == filter)
                .cloned()
        } else {
            names.first().cloned()
        };
        Ok((key, names, selected_enzyme))
    }

    pub(super) fn complement_iupac_text(seq: &str) -> String {
        seq.as_bytes()
            .iter()
            .map(|b| {
                Self::iupac_letter_complement(*b)
                    .unwrap_or(b'N')
                    .to_ascii_uppercase() as char
            })
            .collect()
    }

    pub(super) fn build_restriction_site_expert_view(
        &self,
        seq_id: &str,
        cut_pos_1based: usize,
        enzyme: Option<&str>,
        recognition_start_1based: Option<usize>,
        recognition_end_1based: Option<usize>,
    ) -> Result<RestrictionSiteExpertView, EngineError> {
        let (key, enzyme_names, selected_enzyme) = self.resolve_restriction_site_target(
            seq_id,
            cut_pos_1based,
            enzyme,
            recognition_start_1based,
            recognition_end_1based,
        )?;
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let start_1based = key.from().max(0) as usize + 1;
        let end_1based = key.to().max(0) as usize;
        let start0 = start_1based.saturating_sub(1);
        let end0 = end_1based;
        let mut site_sequence = dna
            .get_range_safe(start0..end0)
            .map(|bytes| String::from_utf8_lossy(&bytes).to_ascii_uppercase())
            .unwrap_or_default();
        let selected_enzyme_def = selected_enzyme.as_ref().and_then(|name| {
            dna.restriction_enzymes()
                .iter()
                .find(|enzyme| enzyme.name.eq_ignore_ascii_case(name))
        });
        let recognition_iupac = selected_enzyme_def.map(|enzyme| enzyme.sequence.clone());
        let enzyme_cut_offset_0based = selected_enzyme_def.map(|enzyme| enzyme.cut);
        let overlap_bp = selected_enzyme_def.map(|enzyme| enzyme.overlap);
        let enzyme_note = selected_enzyme_def.and_then(|enzyme| {
            enzyme
                .note
                .as_deref()
                .map(str::trim)
                .filter(|note| !note.is_empty())
                .map(str::to_string)
        });
        let rebase_url = selected_enzyme_def
            .map(|enzyme| format!("https://rebase.neb.com/rebase/enz/{}.html", enzyme.name));
        if site_sequence.is_empty() {
            if let Some(iupac) = &recognition_iupac {
                site_sequence = iupac.to_ascii_uppercase();
            }
        }
        let site_sequence_complement = Self::complement_iupac_text(&site_sequence);
        let max_cut = site_sequence.chars().count();
        let cut_index_0based = key.cut_size().max(0) as usize;
        let cut_index_0based = cut_index_0based.min(max_cut);

        Ok(RestrictionSiteExpertView {
            seq_id: seq_id.to_string(),
            cut_pos_1based,
            recognition_start_1based: start_1based,
            recognition_end_1based: end_1based,
            cut_index_0based,
            number_of_cuts_for_enzyme: key.number_of_cuts(),
            selected_enzyme,
            enzyme_names,
            recognition_iupac,
            site_sequence,
            site_sequence_complement,
            enzyme_cut_offset_0based,
            overlap_bp,
            enzyme_note,
            rebase_url,
            instruction: RESTRICTION_EXPERT_INSTRUCTION.to_string(),
        })
    }

    pub(super) fn is_mrna_feature(feature: &gb_io::seq::Feature) -> bool {
        let kind = feature.kind.to_string();
        kind.eq_ignore_ascii_case("mRNA") || kind.eq_ignore_ascii_case("transcript")
    }

    pub(super) fn is_splicing_noncoding_transcript_feature(feature: &gb_io::seq::Feature) -> bool {
        let kind = feature.kind.to_string();
        kind.eq_ignore_ascii_case("ncRNA") || kind.eq_ignore_ascii_case("misc_RNA")
    }

    pub(super) fn is_splicing_transcript_feature(feature: &gb_io::seq::Feature) -> bool {
        Self::is_mrna_feature(feature) || Self::is_splicing_noncoding_transcript_feature(feature)
    }

    pub(super) fn is_splicing_seed_feature(feature: &gb_io::seq::Feature) -> bool {
        Self::is_splicing_transcript_feature(feature)
            || Self::is_exon_feature(feature)
            || feature.kind.to_string().eq_ignore_ascii_case("gene")
            || feature.kind.to_string().eq_ignore_ascii_case("cds")
    }

    pub(super) fn is_exon_feature(feature: &gb_io::seq::Feature) -> bool {
        feature.kind.to_string().eq_ignore_ascii_case("exon")
    }

    pub(super) fn splicing_group_label(feature: &gb_io::seq::Feature, fallback: usize) -> String {
        Self::first_nonempty_feature_qualifier(
            feature,
            &[
                "gene",
                "gene_id",
                "locus_tag",
                "standard_name",
                "gene_synonym",
            ],
        )
        .unwrap_or_else(|| format!("mRNA-group-{}", fallback + 1))
    }

    pub(super) fn feature_transcript_id(feature: &gb_io::seq::Feature, fallback: usize) -> String {
        Self::first_nonempty_feature_qualifier(
            feature,
            &[
                "transcript_id",
                "standard_name",
                "product",
                "label",
                "name",
                "gene",
            ],
        )
        .unwrap_or_else(|| format!("transcript-{}", fallback + 1))
    }

    pub(super) fn range_vec_to_splicing(mut ranges: Vec<(usize, usize)>) -> Vec<SplicingRange> {
        ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        ranges
            .into_iter()
            .filter(|(start, end)| end > start)
            .map(|(start, end)| SplicingRange {
                start_1based: start + 1,
                end_1based: end,
            })
            .collect()
    }

    pub(super) fn order_ranges_for_transcript(
        mut ranges: Vec<(usize, usize)>,
        is_reverse: bool,
    ) -> Vec<(usize, usize)> {
        ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        if is_reverse {
            ranges.reverse();
        }
        ranges
    }

    pub(super) fn cds_ranges_to_reference_aa_segments(
        cds_ranges: Vec<(usize, usize)>,
        is_reverse: bool,
        reference_start_aa: Option<usize>,
    ) -> Vec<IsoformArchitectureCdsAaSegment> {
        let ordered = Self::order_ranges_for_transcript(cds_ranges, is_reverse);
        if ordered.is_empty() {
            return vec![];
        }
        let reference_start = reference_start_aa.unwrap_or(1).max(1);
        let mut cumulative_nt = 0usize;
        let mut out = vec![];
        for (start_0based, end_0based_exclusive) in ordered {
            if end_0based_exclusive <= start_0based {
                continue;
            }
            let segment_nt = end_0based_exclusive - start_0based;
            let local_aa_start = cumulative_nt / 3 + 1;
            let local_aa_end = (cumulative_nt + segment_nt) / 3;
            cumulative_nt = cumulative_nt.saturating_add(segment_nt);
            if local_aa_end < local_aa_start {
                continue;
            }
            let aa_start = reference_start.saturating_add(local_aa_start.saturating_sub(1));
            let aa_end = reference_start.saturating_add(local_aa_end.saturating_sub(1));
            out.push(IsoformArchitectureCdsAaSegment {
                genomic_start_1based: start_0based.saturating_add(1),
                genomic_end_1based: end_0based_exclusive,
                aa_start,
                aa_end,
            });
        }
        out
    }

    pub(super) fn sequence_slice_upper(dna: &DNAsequence, start: usize, end: usize) -> String {
        if end <= start {
            return String::new();
        }
        dna.get_range_safe(start..end)
            .map(|bytes| String::from_utf8_lossy(&bytes).to_ascii_uppercase())
            .unwrap_or_default()
    }

    pub(super) fn normalize_transcript_probe(raw: &str) -> String {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return String::new();
        }
        let mut out = trimmed.to_ascii_uppercase();
        if let Some((prefix, suffix)) = out.rsplit_once('.')
            && !prefix.is_empty()
            && !suffix.is_empty()
            && suffix.chars().all(|ch| ch.is_ascii_digit())
        {
            out = prefix.to_string();
        }
        out
    }

    pub(super) fn feature_transcript_match_keys(
        feature: &gb_io::seq::Feature,
        fallback: usize,
    ) -> Vec<String> {
        let mut keys: BTreeSet<String> = BTreeSet::new();
        for key in [
            "transcript_id",
            "standard_name",
            "label",
            "name",
            "product",
            "protein_id",
        ] {
            if let Some(value) = Self::feature_qualifier_text(feature, key) {
                let normalized = Self::normalize_transcript_probe(&value);
                if !normalized.is_empty() {
                    keys.insert(normalized);
                }
            }
        }
        let fallback_id =
            Self::normalize_transcript_probe(&Self::feature_transcript_id(feature, fallback));
        if !fallback_id.is_empty() {
            keys.insert(fallback_id);
        }
        keys.into_iter().collect()
    }

    pub(super) fn normalize_uniprot_projection_id(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "projection_id cannot be empty".to_string(),
            });
        }
        let mut out = String::with_capacity(trimmed.len());
        for ch in trimmed.chars() {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.' | '@') {
                out.push(ch);
            } else {
                out.push('_');
            }
        }
        let normalized = out.trim_matches('_').to_string();
        if normalized.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "projection_id must contain at least one ASCII letter, digit, '.', '-', '_' or '@'"
                        .to_string(),
            });
        }
        Ok(normalized)
    }

    pub(super) fn load_uniprot_entry_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> UniprotEntryStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<UniprotEntryStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = UNIPROT_ENTRIES_SCHEMA.to_string();
        }
        store
    }

    pub(super) fn load_uniprot_projection_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> UniprotGenomeProjectionStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<UniprotGenomeProjectionStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = UNIPROT_GENOME_PROJECTIONS_SCHEMA.to_string();
        }
        store
    }

    pub(super) fn read_uniprot_entry_store(&self) -> UniprotEntryStore {
        Self::load_uniprot_entry_store_from_metadata(
            self.state.metadata.get(UNIPROT_ENTRIES_METADATA_KEY),
        )
    }

    pub(super) fn write_uniprot_entry_store(
        &mut self,
        mut store: UniprotEntryStore,
    ) -> Result<(), EngineError> {
        if store.entries.is_empty() {
            self.state.metadata.remove(UNIPROT_ENTRIES_METADATA_KEY);
            return Ok(());
        }
        store.schema = UNIPROT_ENTRIES_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize UniProt entry metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(UNIPROT_ENTRIES_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub(super) fn read_uniprot_projection_store(&self) -> UniprotGenomeProjectionStore {
        Self::load_uniprot_projection_store_from_metadata(
            self.state
                .metadata
                .get(UNIPROT_GENOME_PROJECTIONS_METADATA_KEY),
        )
    }

    pub(super) fn write_uniprot_projection_store(
        &mut self,
        mut store: UniprotGenomeProjectionStore,
    ) -> Result<(), EngineError> {
        if store.projections.is_empty() {
            self.state
                .metadata
                .remove(UNIPROT_GENOME_PROJECTIONS_METADATA_KEY);
            return Ok(());
        }
        store.schema = UNIPROT_GENOME_PROJECTIONS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize UniProt projection metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(UNIPROT_GENOME_PROJECTIONS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub(super) fn upsert_uniprot_entry(
        &mut self,
        mut entry: UniprotEntry,
    ) -> Result<(), EngineError> {
        let normalized = normalize_uniprot_entry_id(&entry.entry_id);
        if normalized.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "UniProt entry_id cannot be empty".to_string(),
            });
        }
        entry.entry_id = normalized.clone();
        if entry.schema.trim().is_empty() {
            entry.schema = "gentle.uniprot_entry.v1".to_string();
        }
        entry.imported_at_unix_ms = Self::now_unix_ms();
        let mut aliases = BTreeSet::new();
        for value in std::iter::once(entry.entry_id.clone())
            .chain(entry.aliases.clone().into_iter())
            .chain(std::iter::once(entry.accession.clone()))
            .chain(std::iter::once(entry.primary_id.clone()))
        {
            let normalized_alias = normalize_uniprot_entry_id(&value);
            if !normalized_alias.is_empty() {
                aliases.insert(normalized_alias);
            }
        }
        entry.aliases = aliases.into_iter().collect::<Vec<_>>();

        let mut store = self.read_uniprot_entry_store();
        store.entries.insert(entry.entry_id.clone(), entry);
        self.write_uniprot_entry_store(store)
    }

    pub(super) fn get_uniprot_entry_from_store(
        store: &UniprotEntryStore,
        entry_id: &str,
    ) -> Option<UniprotEntry> {
        let probe = normalize_uniprot_entry_id(entry_id);
        if probe.is_empty() {
            return None;
        }
        if let Some(entry) = store.entries.get(&probe) {
            return Some(entry.clone());
        }
        store
            .entries
            .values()
            .find(|entry| {
                entry
                    .aliases
                    .iter()
                    .any(|alias| normalize_uniprot_entry_id(alias) == probe)
            })
            .cloned()
    }

    pub fn get_uniprot_entry(&self, entry_id: &str) -> Result<UniprotEntry, EngineError> {
        let store = self.read_uniprot_entry_store();
        Self::get_uniprot_entry_from_store(&store, entry_id).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "UniProt entry '{}' not found in project metadata",
                entry_id.trim()
            ),
        })
    }

    pub fn list_uniprot_entries(&self) -> Vec<UniprotEntrySummary> {
        let store = self.read_uniprot_entry_store();
        let mut rows = store
            .entries
            .values()
            .map(|entry| UniprotEntrySummary {
                entry_id: entry.entry_id.clone(),
                accession: entry.accession.clone(),
                primary_id: entry.primary_id.clone(),
                sequence_length: entry.sequence_length,
                feature_count: entry.features.len(),
                ensembl_xref_count: entry.ensembl_xrefs.len(),
                nucleotide_xref_count: entry.nucleotide_xrefs.len(),
                imported_at_unix_ms: entry.imported_at_unix_ms,
                source: entry.source.clone(),
            })
            .collect::<Vec<_>>();
        rows.sort_by(|a, b| a.entry_id.cmp(&b.entry_id));
        rows
    }

    pub fn get_uniprot_genome_projection(
        &self,
        projection_id: &str,
    ) -> Result<UniprotGenomeProjection, EngineError> {
        let projection_id = Self::normalize_uniprot_projection_id(projection_id)?;
        let store = self.read_uniprot_projection_store();
        store
            .projections
            .get(&projection_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "UniProt genome projection '{}' was not found",
                    projection_id
                ),
            })
    }

    pub fn list_uniprot_genome_projections(
        &self,
        seq_filter: Option<&str>,
    ) -> Vec<UniprotGenomeProjectionSummary> {
        let seq_filter = seq_filter
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_string());
        let store = self.read_uniprot_projection_store();
        let mut rows = store
            .projections
            .values()
            .filter(|record| {
                seq_filter
                    .as_deref()
                    .map(|seq_id| record.seq_id == seq_id)
                    .unwrap_or(true)
            })
            .map(|record| UniprotGenomeProjectionSummary {
                projection_id: record.projection_id.clone(),
                entry_id: record.entry_id.clone(),
                seq_id: record.seq_id.clone(),
                created_at_unix_ms: record.created_at_unix_ms,
                transcript_projection_count: record.transcript_projections.len(),
            })
            .collect::<Vec<_>>();
        rows.sort_by(|a, b| a.projection_id.cmp(&b.projection_id));
        rows
    }

    pub(super) fn fetch_uniprot_swiss_prot_text(
        query: &str,
    ) -> Result<(String, String), EngineError> {
        let query = query.trim();
        if query.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "UniProt query cannot be empty".to_string(),
            });
        }
        let query_url = query.replace(' ', "%20");
        let url = format!("https://rest.uniprot.org/uniprotkb/{query_url}.txt");
        let client = reqwest::blocking::Client::builder()
            .timeout(Duration::from_secs(45))
            .build()
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not create UniProt HTTP client: {e}"),
            })?;
        let response = client.get(&url).send().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not fetch UniProt entry '{query}': {e}"),
        })?;
        if !response.status().is_success() {
            let status = response.status();
            let detail = if status == reqwest::StatusCode::NOT_FOUND {
                format!(
                    "UniProt query '{}' was not found (HTTP {}) at '{}'; verify accession/entry ID",
                    query, status, url
                )
            } else {
                format!(
                    "UniProt query '{}' returned HTTP status {} at '{}'",
                    query, status, url
                )
            };
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: detail,
            });
        }
        let text = response.text().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read UniProt response body for '{query}': {e}"),
        })?;
        Ok((url, text))
    }

    pub(super) fn fetch_genbank_accession_text(
        accession: &str,
    ) -> Result<(String, String), EngineError> {
        let accession = validate_genbank_accession(accession).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: e,
        })?;
        let source_url = build_genbank_efetch_url(&accession, "gbwithparts");
        if let Some(path) = source_url.strip_prefix("file://") {
            let text = std::fs::read_to_string(path).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not read GenBank source file '{}' for accession '{}': {e}",
                    path, accession
                ),
            })?;
            return Ok((source_url, text));
        }
        let client = reqwest::blocking::Client::builder()
            .timeout(Duration::from_secs(45))
            .build()
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not create GenBank HTTP client: {e}"),
            })?;
        let response = client.get(&source_url).send().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not fetch GenBank accession '{}' from '{}': {e}",
                accession, source_url
            ),
        })?;
        if !response.status().is_success() {
            let status = response.status();
            let detail = if status == reqwest::StatusCode::NOT_FOUND {
                format!(
                    "GenBank accession '{}' was not found (HTTP {}) at '{}'",
                    accession, status, source_url
                )
            } else {
                format!(
                    "GenBank accession '{}' returned HTTP status {} at '{}'",
                    accession, status, source_url
                )
            };
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: detail,
            });
        }
        let text = response.text().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not read GenBank response body for accession '{}': {e}",
                accession
            ),
        })?;
        Ok((source_url, text))
    }

    pub(super) fn parse_uniprot_entry_text(
        text: &str,
        source: &str,
        source_query: Option<&str>,
        entry_id_override: Option<&str>,
    ) -> Result<UniprotEntry, EngineError> {
        parse_swiss_prot_text(text, source, source_query, entry_id_override).map_err(|e| {
            EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not parse SWISS-PROT text: {e}"),
            }
        })
    }

    pub(super) fn import_uniprot_entry_sequence(
        &mut self,
        _result: &mut OpResult,
        entry_id: &str,
        _output_id: Option<&str>,
    ) -> Result<SeqId, EngineError> {
        let _entry = self.get_uniprot_entry(entry_id)?;
        Err(EngineError {
            code: ErrorCode::Unsupported,
            message: format!(
                "ImportUniprotEntrySequence is currently disabled: UniProt protein sequences are metadata/projection-only in this release (entry_id='{}').",
                entry_id
            ),
        })
    }

    pub(super) fn aa_interval_to_genomic_segments(
        aa_start: usize,
        aa_end: usize,
        aa_segments: &[UniprotAaGenomicSegment],
        is_reverse: bool,
    ) -> Vec<UniprotAaGenomicSegment> {
        let mut out = vec![];
        for segment in aa_segments {
            if segment.aa_end < aa_start || segment.aa_start > aa_end {
                continue;
            }
            let overlap_start = aa_start.max(segment.aa_start);
            let overlap_end = aa_end.min(segment.aa_end);
            if overlap_end < overlap_start {
                continue;
            }
            let seg_nt_start = segment.genomic_start_1based.min(segment.genomic_end_1based);
            let seg_nt_end = segment.genomic_start_1based.max(segment.genomic_end_1based);
            let seg_nt_len = seg_nt_end.saturating_sub(seg_nt_start).saturating_add(1);
            let seg_aa_len = segment
                .aa_end
                .saturating_sub(segment.aa_start)
                .saturating_add(1);
            if seg_nt_len == 0 || seg_aa_len == 0 {
                continue;
            }
            let nt_start_offset = overlap_start
                .saturating_sub(segment.aa_start)
                .saturating_mul(3);
            let nt_end_offset_exclusive = overlap_end
                .saturating_sub(segment.aa_start)
                .saturating_add(1)
                .saturating_mul(3);
            let nt_start_offset = nt_start_offset.min(seg_nt_len.saturating_sub(1));
            let nt_end_offset_exclusive = nt_end_offset_exclusive.min(seg_nt_len);
            if nt_end_offset_exclusive <= nt_start_offset {
                continue;
            }
            let (genomic_start_1based, genomic_end_1based) = if is_reverse {
                let end = seg_nt_end.saturating_sub(nt_start_offset);
                let start = seg_nt_end.saturating_sub(nt_end_offset_exclusive.saturating_sub(1));
                (start.min(end), start.max(end))
            } else {
                let start = seg_nt_start.saturating_add(nt_start_offset);
                let end = seg_nt_start.saturating_add(nt_end_offset_exclusive.saturating_sub(1));
                (start.min(end), start.max(end))
            };
            out.push(UniprotAaGenomicSegment {
                aa_start: overlap_start,
                aa_end: overlap_end,
                genomic_start_1based,
                genomic_end_1based,
                strand: segment.strand.clone(),
            });
        }
        out
    }

    pub fn project_uniprot_to_genome(
        &self,
        seq_id: &str,
        entry_id: &str,
        projection_id: Option<&str>,
        transcript_id_filter: Option<&str>,
    ) -> Result<UniprotGenomeProjection, EngineError> {
        let seq = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", seq_id),
            })?;
        let entry = self.get_uniprot_entry(entry_id)?;
        let mut warnings: Vec<String> = vec![];
        let mut transcripts_requested: Vec<String> = vec![];
        if let Some(filter) = transcript_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            transcripts_requested.push(filter.to_string());
        } else {
            for xref in &entry.ensembl_xrefs {
                if let Some(transcript_id) = xref
                    .transcript_id
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    transcripts_requested.push(transcript_id.to_string());
                }
            }
        }

        let mut transcript_feature_keys: HashMap<String, usize> = HashMap::new();
        let features = seq.features();
        for (idx, feature) in features.iter().enumerate() {
            if !Self::is_mrna_feature(feature) {
                continue;
            }
            for key in Self::feature_transcript_match_keys(feature, idx) {
                transcript_feature_keys.entry(key).or_insert(idx);
            }
        }

        if transcripts_requested.is_empty() {
            warnings.push(
                "UniProt entry has no transcript xref; falling back to all mRNA features in sequence"
                    .to_string(),
            );
            for (idx, feature) in features.iter().enumerate() {
                if !Self::is_mrna_feature(feature) {
                    continue;
                }
                transcripts_requested.push(Self::feature_transcript_id(feature, idx));
            }
        }

        let mut transcript_projections: Vec<UniprotTranscriptProjection> = vec![];
        let mut seen_transcripts = BTreeSet::new();
        for transcript_probe in transcripts_requested {
            let normalized_probe = Self::normalize_transcript_probe(&transcript_probe);
            if normalized_probe.is_empty() || !seen_transcripts.insert(normalized_probe.clone()) {
                continue;
            }
            let Some(feature_index) = transcript_feature_keys.get(&normalized_probe).copied()
            else {
                warnings.push(format!(
                    "Transcript '{}' from UniProt entry '{}' was not found in sequence '{}'",
                    transcript_probe, entry.entry_id, seq_id
                ));
                continue;
            };
            let feature = &features[feature_index];
            let transcript_id = Self::feature_transcript_id(feature, feature_index);
            let is_reverse = feature_is_reverse(feature);
            let mut transcript_warnings: Vec<String> = vec![];
            let mut cds_ranges =
                Self::feature_qualifier_ranges_0based(feature, "cds_ranges_1based");
            if cds_ranges.is_empty() {
                let mut fallback = vec![];
                collect_location_ranges_usize(&feature.location, &mut fallback);
                fallback.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
                fallback.retain(|(start, end)| end > start);
                cds_ranges = fallback;
                transcript_warnings.push(format!(
                    "Transcript '{}' has no cds_ranges_1based qualifier; fell back to feature exon spans",
                    transcript_id
                ));
            }
            let aa_reference_segments =
                Self::cds_ranges_to_reference_aa_segments(cds_ranges.clone(), is_reverse, Some(1))
                    .into_iter()
                    .map(|segment| UniprotAaGenomicSegment {
                        aa_start: segment.aa_start,
                        aa_end: segment.aa_end,
                        genomic_start_1based: segment.genomic_start_1based,
                        genomic_end_1based: segment.genomic_end_1based,
                        strand: if is_reverse {
                            "-".to_string()
                        } else {
                            "+".to_string()
                        },
                    })
                    .collect::<Vec<_>>();
            if aa_reference_segments.is_empty() {
                transcript_warnings.push(format!(
                    "Transcript '{}' has no mappable CDS amino-acid segments",
                    transcript_id
                ));
            }

            let mut feature_projections: Vec<UniprotFeatureProjection> = vec![];
            for feature in &entry.features {
                let (Some(aa_start), Some(aa_end)) =
                    (feature.interval.start_aa, feature.interval.end_aa)
                else {
                    continue;
                };
                if aa_end < aa_start {
                    continue;
                }
                let genomic_segments = Self::aa_interval_to_genomic_segments(
                    aa_start,
                    aa_end,
                    &aa_reference_segments,
                    is_reverse,
                );
                if genomic_segments.is_empty() {
                    continue;
                }
                feature_projections.push(UniprotFeatureProjection {
                    feature_key: feature.key.clone(),
                    feature_note: feature.note.clone(),
                    aa_start,
                    aa_end,
                    genomic_segments,
                });
            }
            transcript_projections.push(UniprotTranscriptProjection {
                transcript_id,
                transcript_feature_id: Some(feature_index),
                strand: if is_reverse {
                    "-".to_string()
                } else {
                    "+".to_string()
                },
                aa_segments: aa_reference_segments,
                feature_projections,
                warnings: transcript_warnings,
            });
        }

        let projection_id = projection_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string())
            .unwrap_or_else(|| format!("{}@{}", entry.entry_id, seq_id));
        let projection_id = Self::normalize_uniprot_projection_id(&projection_id)?;
        Ok(UniprotGenomeProjection {
            schema: UNIPROT_GENOME_PROJECTION_SCHEMA.to_string(),
            projection_id,
            entry_id: entry.entry_id.clone(),
            seq_id: seq_id.to_string(),
            created_at_unix_ms: Self::now_unix_ms(),
            transcript_id_filter: transcript_id_filter
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            transcript_projections,
            warnings,
        })
    }

    pub(super) fn upsert_uniprot_projection(
        &mut self,
        projection: UniprotGenomeProjection,
    ) -> Result<(), EngineError> {
        let mut store = self.read_uniprot_projection_store();
        store
            .projections
            .insert(projection.projection_id.clone(), projection);
        self.write_uniprot_projection_store(store)
    }

    pub(super) fn load_isoform_panel_resource(
        path: &str,
        panel_id_override: Option<&str>,
    ) -> Result<IsoformPanelResource, EngineError> {
        let text = std::fs::read_to_string(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read isoform panel file '{path}': {e}"),
        })?;
        let mut resource =
            serde_json::from_str::<IsoformPanelResource>(&text).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not parse isoform panel JSON in '{path}': {e}"),
            })?;
        if resource.schema.trim().is_empty() {
            resource.schema = ISOFORM_PANEL_RESOURCE_SCHEMA.to_string();
        }
        if resource.schema != ISOFORM_PANEL_RESOURCE_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported isoform-panel schema '{}' in '{}'; expected '{}'",
                    resource.schema, path, ISOFORM_PANEL_RESOURCE_SCHEMA
                ),
            });
        }
        if let Some(override_id) = panel_id_override.map(str::trim).filter(|v| !v.is_empty()) {
            resource.panel_id = override_id.to_string();
        }
        resource.panel_id = resource.panel_id.trim().to_string();
        if resource.panel_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Isoform panel file '{}' has empty panel_id", path),
            });
        }
        resource.gene_symbol = resource.gene_symbol.trim().to_string();
        if resource.gene_symbol.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Isoform panel file '{}' has empty gene_symbol", path),
            });
        }
        if resource.isoforms.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Isoform panel file '{}' has no isoforms", path),
            });
        }
        for (idx, isoform) in resource.isoforms.iter_mut().enumerate() {
            isoform.isoform_id = isoform.isoform_id.trim().to_string();
            if isoform.isoform_id.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Isoform panel '{}' contains isoform with empty isoform_id at row {}",
                        resource.panel_id,
                        idx + 1
                    ),
                });
            }
            isoform.label = isoform
                .label
                .take()
                .map(|v| v.trim().to_string())
                .filter(|v| !v.is_empty());
            let mut transcript_ids: Vec<String> = vec![];
            let mut seen_transcripts: BTreeSet<String> = BTreeSet::new();
            for entry in &isoform.transcript_ids {
                let normalized = entry.trim().to_string();
                if normalized.is_empty() {
                    continue;
                }
                let key = normalized.to_ascii_uppercase();
                if seen_transcripts.insert(key) {
                    transcript_ids.push(normalized);
                }
            }
            isoform.transcript_ids = transcript_ids;
            for (domain_idx, domain) in isoform.domains.iter_mut().enumerate() {
                domain.name = domain.name.trim().to_string();
                if domain.name.is_empty() {
                    domain.name = format!("domain_{}", domain_idx + 1);
                }
                if domain.start_aa == 0 || domain.end_aa == 0 || domain.end_aa < domain.start_aa {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Isoform panel '{}' has invalid domain range for isoform '{}' domain '{}': {}..{}",
                            resource.panel_id,
                            isoform.isoform_id,
                            domain.name,
                            domain.start_aa,
                            domain.end_aa
                        ),
                    });
                }
                domain.color_hex = domain
                    .color_hex
                    .take()
                    .map(|v| v.trim().to_string())
                    .filter(|v| !v.is_empty());
            }
            isoform.domains.sort_by(|left, right| {
                left.start_aa
                    .cmp(&right.start_aa)
                    .then(left.end_aa.cmp(&right.end_aa))
            });
        }
        Ok(resource)
    }

    pub(super) fn isoform_validation_issue(
        code: &str,
        message: String,
        isoform_id: Option<&str>,
        transcript_probe: Option<&str>,
        domain_name: Option<&str>,
    ) -> IsoformPanelValidationIssue {
        IsoformPanelValidationIssue {
            severity: "warning".to_string(),
            code: code.to_string(),
            message,
            isoform_id: isoform_id.map(|v| v.to_string()),
            transcript_probe: transcript_probe.map(|v| v.to_string()),
            domain_name: domain_name.map(|v| v.to_string()),
        }
    }

    pub fn validate_isoform_panel_resource(
        path: &str,
        panel_id_override: Option<&str>,
    ) -> Result<IsoformPanelValidationReport, EngineError> {
        let resource = Self::load_isoform_panel_resource(path, panel_id_override)?;
        let color_hex_re = Regex::new(r"^#[0-9A-Fa-f]{6}$").map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not compile isoform color-hex regex: {e}"),
        })?;

        let mut issues: Vec<IsoformPanelValidationIssue> = vec![];
        let mut isoforms: Vec<IsoformPanelValidationIsoformSummary> = vec![];
        let mut transcript_probe_count = 0usize;
        let mut unique_transcript_probes: BTreeSet<String> = BTreeSet::new();
        let mut domain_count = 0usize;
        let mut isoform_id_buckets: HashMap<String, Vec<String>> = HashMap::new();
        let mut transcript_probe_buckets: HashMap<String, BTreeSet<String>> = HashMap::new();

        for isoform in &resource.isoforms {
            let isoform_key = isoform.isoform_id.to_ascii_uppercase();
            isoform_id_buckets
                .entry(isoform_key)
                .or_default()
                .push(isoform.isoform_id.clone());

            if isoform.transcript_ids.is_empty() {
                issues.push(Self::isoform_validation_issue(
                    "missing_transcript_probe",
                    format!(
                        "Isoform '{}' has no transcript probes; mapping can only rely on label heuristics",
                        isoform.isoform_id
                    ),
                    Some(&isoform.isoform_id),
                    None,
                    None,
                ));
            }
            for transcript_probe in &isoform.transcript_ids {
                let normalized = Self::normalize_transcript_probe(transcript_probe);
                if normalized.is_empty() {
                    continue;
                }
                transcript_probe_count += 1;
                unique_transcript_probes.insert(normalized.clone());
                transcript_probe_buckets
                    .entry(normalized)
                    .or_default()
                    .insert(isoform.isoform_id.clone());
            }

            if isoform.domains.is_empty() {
                issues.push(Self::isoform_validation_issue(
                    "missing_domains",
                    format!(
                        "Isoform '{}' has no domains; protein lane may render without functional landmarks",
                        isoform.isoform_id
                    ),
                    Some(&isoform.isoform_id),
                    None,
                    None,
                ));
            }

            domain_count += isoform.domains.len();
            let mut sorted_domains = isoform.domains.clone();
            sorted_domains.sort_by(|left, right| {
                left.start_aa
                    .cmp(&right.start_aa)
                    .then(left.end_aa.cmp(&right.end_aa))
                    .then(left.name.cmp(&right.name))
            });

            for domain in &sorted_domains {
                if let Some(color_hex) = domain.color_hex.as_deref()
                    && !color_hex_re.is_match(color_hex)
                {
                    issues.push(Self::isoform_validation_issue(
                        "invalid_color_hex",
                        format!(
                            "Isoform '{}' domain '{}' uses non-hex color '{}'; expected '#RRGGBB'",
                            isoform.isoform_id, domain.name, color_hex
                        ),
                        Some(&isoform.isoform_id),
                        None,
                        Some(&domain.name),
                    ));
                }
            }

            for pair in sorted_domains.windows(2) {
                if pair[1].start_aa <= pair[0].end_aa {
                    issues.push(Self::isoform_validation_issue(
                        "overlapping_domains",
                        format!(
                            "Isoform '{}' has overlapping domains '{}' ({}..{}) and '{}' ({}..{})",
                            isoform.isoform_id,
                            pair[0].name,
                            pair[0].start_aa,
                            pair[0].end_aa,
                            pair[1].name,
                            pair[1].start_aa,
                            pair[1].end_aa
                        ),
                        Some(&isoform.isoform_id),
                        None,
                        Some(&pair[1].name),
                    ));
                }
            }

            let max_domain_end_aa = sorted_domains.iter().map(|domain| domain.end_aa).max();
            if let (Some(reference_start), Some(reference_end)) =
                (isoform.reference_start_aa, isoform.reference_end_aa)
                && reference_end < reference_start
            {
                issues.push(Self::isoform_validation_issue(
                    "invalid_reference_range",
                    format!(
                        "Isoform '{}' reference range is invalid (start={} > end={})",
                        isoform.isoform_id, reference_start, reference_end
                    ),
                    Some(&isoform.isoform_id),
                    None,
                    None,
                ));
            }
            if let Some(reference_end) = isoform.reference_end_aa {
                if let Some(max_domain_end) = max_domain_end_aa
                    && reference_end < max_domain_end
                {
                    issues.push(Self::isoform_validation_issue(
                        "reference_end_below_domain_end",
                        format!(
                            "Isoform '{}' reference_end_aa={} is smaller than max domain end {}",
                            isoform.isoform_id, reference_end, max_domain_end
                        ),
                        Some(&isoform.isoform_id),
                        None,
                        None,
                    ));
                }
            } else if let (Some(expected), Some(max_domain_end)) =
                (isoform.expected_length_aa, max_domain_end_aa)
                && expected < max_domain_end
            {
                issues.push(Self::isoform_validation_issue(
                    "expected_length_below_domain_end",
                    format!(
                        "Isoform '{}' expected_length_aa={} is smaller than max domain end {}",
                        isoform.isoform_id, expected, max_domain_end
                    ),
                    Some(&isoform.isoform_id),
                    None,
                    None,
                ));
            }

            isoforms.push(IsoformPanelValidationIsoformSummary {
                isoform_id: isoform.isoform_id.clone(),
                label: isoform
                    .label
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or(isoform.isoform_id.as_str())
                    .to_string(),
                transcript_probe_count: isoform.transcript_ids.len(),
                domain_count: isoform.domains.len(),
                expected_length_aa: isoform.expected_length_aa,
                max_domain_end_aa,
            });
        }

        for isoform_ids in isoform_id_buckets.values() {
            if isoform_ids.len() < 2 {
                continue;
            }
            let mut ids = isoform_ids.clone();
            ids.sort();
            ids.dedup();
            issues.push(Self::isoform_validation_issue(
                "duplicate_isoform_id",
                format!(
                    "Panel '{}' contains duplicate isoform_id entries: {}",
                    resource.panel_id,
                    ids.join(", ")
                ),
                ids.first().map(String::as_str),
                None,
                None,
            ));
        }

        for (normalized_probe, isoform_ids) in &transcript_probe_buckets {
            if isoform_ids.len() < 2 {
                continue;
            }
            let ids = isoform_ids.iter().cloned().collect::<Vec<_>>();
            issues.push(Self::isoform_validation_issue(
                "shared_transcript_probe",
                format!(
                    "Transcript probe '{}' is shared by multiple isoforms: {}",
                    normalized_probe,
                    ids.join(", ")
                ),
                ids.first().map(String::as_str),
                Some(normalized_probe.as_str()),
                None,
            ));
        }

        issues.sort_by(|left, right| {
            left.code
                .cmp(&right.code)
                .then(left.isoform_id.cmp(&right.isoform_id))
                .then(left.transcript_probe.cmp(&right.transcript_probe))
                .then(left.domain_name.cmp(&right.domain_name))
                .then(left.message.cmp(&right.message))
        });

        let issue_count = issues.len();
        let status = if issue_count == 0 {
            "ok".to_string()
        } else {
            "warning".to_string()
        };

        Ok(IsoformPanelValidationReport {
            schema: ISOFORM_PANEL_VALIDATION_REPORT_SCHEMA.to_string(),
            path: path.to_string(),
            panel_id: resource.panel_id,
            gene_symbol: resource.gene_symbol,
            assembly: resource.assembly,
            isoform_count: isoforms.len(),
            transcript_probe_count,
            unique_transcript_probe_count: unique_transcript_probes.len(),
            domain_count,
            issue_count,
            status,
            isoforms,
            issues,
        })
    }

    pub(super) fn build_isoform_architecture_expert_view_from_resource(
        &self,
        seq_id: &str,
        panel_id: &str,
        resource: &IsoformPanelResource,
        strict: bool,
    ) -> Result<IsoformArchitectureExpertView, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let features = dna.features();

        #[derive(Clone)]
        struct TranscriptProjection {
            feature_id: usize,
            transcript_id: String,
            is_reverse: bool,
            full_exon_ranges: Vec<(usize, usize)>,
            geometry_ranges: Vec<(usize, usize)>,
            cds_ranges: Vec<(usize, usize)>,
            intron_ranges: Vec<(usize, usize)>,
            used_cds_geometry: bool,
        }

        let mut projections: Vec<TranscriptProjection> = vec![];
        let mut projection_key_index: HashMap<String, Vec<usize>> = HashMap::new();
        let use_cds_geometry =
            resource.transcript_geometry_mode == IsoformTranscriptGeometryMode::Cds;

        let gene_probe = resource.gene_symbol.trim().to_ascii_uppercase();
        let mut has_gene_probe_match = gene_probe.is_empty();
        for (feature_id, feature) in features.iter().enumerate() {
            if !has_gene_probe_match
                && let Some(gene_name) = Self::first_nonempty_feature_qualifier(
                    feature,
                    &[
                        "gene",
                        "gene_id",
                        "locus_tag",
                        "standard_name",
                        "name",
                        "product",
                    ],
                )
            {
                has_gene_probe_match = gene_name.to_ascii_uppercase().contains(&gene_probe);
            }

            if !Self::is_mrna_feature(feature) {
                continue;
            }

            let mut exon_ranges = vec![];
            collect_location_ranges_usize(&feature.location, &mut exon_ranges);
            if exon_ranges.is_empty() {
                let (from, to) = feature.location.find_bounds().map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Could not parse transcript range: {e}"),
                })?;
                if from >= 0 && to >= 0 {
                    exon_ranges.push((from as usize, to as usize));
                }
            }
            exon_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            exon_ranges.retain(|(start, end)| end > start);
            if exon_ranges.is_empty() {
                continue;
            }
            let cds_ranges = Self::feature_qualifier_ranges_0based(feature, "cds_ranges_1based");
            let full_exon_ranges = exon_ranges.clone();
            let (geometry_ranges, used_cds_geometry) = if use_cds_geometry && !cds_ranges.is_empty()
            {
                (cds_ranges.clone(), true)
            } else {
                (exon_ranges, false)
            };
            let mut intron_ranges = vec![];
            for pair in geometry_ranges.windows(2) {
                if pair[1].0 > pair[0].1 {
                    intron_ranges.push((pair[0].1, pair[1].0));
                }
            }
            let projection_idx = projections.len();
            let transcript_id = Self::feature_transcript_id(feature, feature_id);
            projections.push(TranscriptProjection {
                feature_id,
                transcript_id: transcript_id.clone(),
                is_reverse: feature_is_reverse(feature),
                full_exon_ranges,
                geometry_ranges,
                cds_ranges,
                intron_ranges,
                used_cds_geometry,
            });
            for key in Self::feature_transcript_match_keys(feature, feature_id) {
                projection_key_index
                    .entry(key)
                    .or_default()
                    .push(projection_idx);
            }
        }

        let mut warnings: Vec<String> = vec![];
        if projections.is_empty() {
            let message = format!(
                "No mRNA/transcript features available in sequence '{}' for isoform panel '{}'",
                seq_id, panel_id
            );
            if strict {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message,
                });
            }
            warnings.push(message);
        }

        if !has_gene_probe_match {
            let message = format!(
                "No feature-gene qualifier in '{}' matched panel gene_symbol '{}'",
                seq_id, resource.gene_symbol
            );
            if strict {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message,
                });
            }
            warnings.push(message);
        }

        let mut mapped_ranges: Vec<(usize, usize)> = vec![];
        let mut transcript_lanes: Vec<IsoformArchitectureTranscriptLane> = vec![];
        let mut protein_lanes: Vec<IsoformArchitectureProteinLane> = vec![];

        for isoform in &resource.isoforms {
            let label = isoform
                .label
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .unwrap_or(isoform.isoform_id.as_str())
                .to_string();
            let probes: Vec<String> = isoform
                .transcript_ids
                .iter()
                .map(|entry| Self::normalize_transcript_probe(entry))
                .filter(|entry| !entry.is_empty())
                .collect();
            let mut matched_projection: Option<&TranscriptProjection> = None;
            for probe in &probes {
                if let Some(indices) = projection_key_index.get(probe)
                    && let Some(first) = indices.first()
                {
                    matched_projection = projections.get(*first);
                    if matched_projection.is_some() {
                        break;
                    }
                }
            }

            let (
                mapped,
                transcript_id,
                transcript_feature_id,
                strand,
                transcript_exons,
                exons,
                introns,
                note,
                cds_to_protein_segments,
            ) = if let Some(projection) = matched_projection {
                mapped_ranges.extend(projection.geometry_ranges.iter().copied());
                let geometry_note = (use_cds_geometry && !projection.used_cds_geometry)
                    .then_some("CDS geometry unavailable; fell back to exon geometry".to_string());
                let cds_to_protein_segments = Self::cds_ranges_to_reference_aa_segments(
                    projection.cds_ranges.clone(),
                    projection.is_reverse,
                    isoform.reference_start_aa,
                );
                (
                    true,
                    Some(projection.transcript_id.clone()),
                    Some(projection.feature_id),
                    if projection.is_reverse {
                        "-".to_string()
                    } else {
                        "+".to_string()
                    },
                    Self::range_vec_to_splicing(projection.full_exon_ranges.clone()),
                    Self::range_vec_to_splicing(projection.geometry_ranges.clone()),
                    Self::range_vec_to_splicing(projection.intron_ranges.clone()),
                    geometry_note,
                    cds_to_protein_segments,
                )
            } else {
                let probe_text = if isoform.transcript_ids.is_empty() {
                    "(none)".to_string()
                } else {
                    isoform.transcript_ids.join(", ")
                };
                let message = format!(
                    "Isoform '{}' has no mapped transcript in '{}' (probes: {})",
                    isoform.isoform_id, seq_id, probe_text
                );
                if strict {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message,
                    });
                }
                warnings.push(message.clone());
                (
                    false,
                    isoform.transcript_ids.first().cloned(),
                    None,
                    "?".to_string(),
                    vec![],
                    vec![],
                    vec![],
                    Some("No transcript mapping found".to_string()),
                    vec![],
                )
            };

            transcript_lanes.push(IsoformArchitectureTranscriptLane {
                isoform_id: isoform.isoform_id.clone(),
                label: label.clone(),
                transcript_id: transcript_id.clone(),
                transcript_feature_id,
                strand,
                transcript_exons,
                exons,
                introns,
                mapped,
                transactivation_class: isoform.transactivation_class.clone(),
                cds_to_protein_segments,
                note,
            });

            let mut domains = isoform
                .domains
                .iter()
                .map(|domain| IsoformArchitectureProteinDomain {
                    name: domain.name.clone(),
                    start_aa: domain.start_aa,
                    end_aa: domain.end_aa,
                    color_hex: domain.color_hex.clone(),
                })
                .collect::<Vec<_>>();
            domains.sort_by(|left, right| {
                left.start_aa
                    .cmp(&right.start_aa)
                    .then(left.end_aa.cmp(&right.end_aa))
            });
            let expected_length_aa = isoform.expected_length_aa.or_else(|| {
                domains
                    .iter()
                    .map(|domain| domain.end_aa)
                    .max()
                    .filter(|max_end| *max_end > 0)
            });
            protein_lanes.push(IsoformArchitectureProteinLane {
                isoform_id: isoform.isoform_id.clone(),
                label,
                transcript_id,
                expected_length_aa,
                reference_start_aa: isoform.reference_start_aa,
                reference_end_aa: isoform.reference_end_aa,
                domains,
                transactivation_class: isoform.transactivation_class.clone(),
            });
        }

        let (region_start_1based, region_end_1based) = if mapped_ranges.is_empty() {
            (1, dna.len().max(1))
        } else {
            let start0 = mapped_ranges
                .iter()
                .map(|(start, _)| *start)
                .min()
                .unwrap_or(0);
            let end0 = mapped_ranges
                .iter()
                .map(|(_, end)| *end)
                .max()
                .unwrap_or(start0 + 1);
            (start0 + 1, end0.max(start0 + 1))
        };

        Ok(IsoformArchitectureExpertView {
            seq_id: seq_id.to_string(),
            panel_id: panel_id.to_string(),
            gene_symbol: resource.gene_symbol.clone(),
            transcript_geometry_mode: resource.transcript_geometry_mode.as_str().to_string(),
            panel_source: resource.source.clone(),
            region_start_1based,
            region_end_1based,
            instruction: ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION.to_string(),
            transcript_lanes,
            protein_lanes,
            warnings,
        })
    }

    pub(super) fn build_isoform_architecture_expert_view(
        &self,
        seq_id: &str,
        panel_id: &str,
    ) -> Result<IsoformArchitectureExpertView, EngineError> {
        let record = self.get_isoform_panel_record(seq_id, panel_id)?;
        self.build_isoform_architecture_expert_view_from_resource(
            seq_id,
            &record.panel_id,
            &record.resource,
            false,
        )
    }

    pub(super) fn splice_boundary_markers_for_introns(
        dna: &DNAsequence,
        transcript_feature_id: usize,
        transcript_id: &str,
        is_reverse: bool,
        introns: &[(usize, usize)],
    ) -> Vec<SplicingBoundaryMarker> {
        let mut out = Vec::new();
        for (start, end) in introns {
            if *end <= *start || end - start < 2 {
                continue;
            }
            if is_reverse {
                let donor_raw = Self::sequence_slice_upper(dna, end.saturating_sub(2), *end);
                let acceptor_raw = Self::sequence_slice_upper(dna, *start, (start + 2).min(*end));
                let donor_rc = Self::reverse_complement_bytes(donor_raw.as_bytes());
                let acceptor_rc = Self::reverse_complement_bytes(acceptor_raw.as_bytes());
                let donor_motif = String::from_utf8_lossy(&donor_rc).to_ascii_uppercase();
                let acceptor_motif = String::from_utf8_lossy(&acceptor_rc).to_ascii_uppercase();
                out.push(SplicingBoundaryMarker {
                    transcript_feature_id,
                    transcript_id: transcript_id.to_string(),
                    side: "donor".to_string(),
                    position_1based: *end,
                    canonical: donor_motif == "GT",
                    motif_2bp: donor_motif,
                });
                out.push(SplicingBoundaryMarker {
                    transcript_feature_id,
                    transcript_id: transcript_id.to_string(),
                    side: "acceptor".to_string(),
                    position_1based: start + 1,
                    canonical: acceptor_motif == "AG",
                    motif_2bp: acceptor_motif,
                });
            } else {
                let donor_motif = Self::sequence_slice_upper(dna, *start, (start + 2).min(*end));
                let acceptor_motif = Self::sequence_slice_upper(dna, end.saturating_sub(2), *end);
                out.push(SplicingBoundaryMarker {
                    transcript_feature_id,
                    transcript_id: transcript_id.to_string(),
                    side: "donor".to_string(),
                    position_1based: start + 1,
                    canonical: donor_motif == "GT",
                    motif_2bp: donor_motif,
                });
                out.push(SplicingBoundaryMarker {
                    transcript_feature_id,
                    transcript_id: transcript_id.to_string(),
                    side: "acceptor".to_string(),
                    position_1based: *end,
                    canonical: acceptor_motif == "AG",
                    motif_2bp: acceptor_motif,
                });
            }
        }
        out
    }

    pub(super) fn range_intersection_0based(
        left: (usize, usize),
        right: (usize, usize),
    ) -> Option<(usize, usize)> {
        let start = left.0.max(right.0);
        let end = left.1.min(right.1);
        (end > start).then_some((start, end))
    }

    pub(super) fn exon_cds_phases_for_transcript(
        feature: &gb_io::seq::Feature,
        exon_ranges: &[(usize, usize)],
        is_reverse: bool,
    ) -> Vec<SplicingExonCdsPhase> {
        let mut phases = exon_ranges
            .iter()
            .map(|(start, end)| SplicingExonCdsPhase {
                start_1based: start + 1,
                end_1based: *end,
                left_cds_phase: None,
                right_cds_phase: None,
            })
            .collect::<Vec<_>>();
        if exon_ranges.is_empty() {
            return phases;
        }
        let mut cds_ranges = Self::feature_qualifier_ranges_0based(feature, "cds_ranges_1based");
        if cds_ranges.is_empty() {
            return phases;
        }
        cds_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        cds_ranges.dedup();

        let mut consumed_cds_bp = 0usize;
        let exon_indices = if is_reverse {
            (0..exon_ranges.len()).rev().collect::<Vec<_>>()
        } else {
            (0..exon_ranges.len()).collect::<Vec<_>>()
        };

        for exon_idx in exon_indices {
            let exon = exon_ranges[exon_idx];
            let mut coding_segments = cds_ranges
                .iter()
                .filter_map(|cds| Self::range_intersection_0based(exon, *cds))
                .collect::<Vec<_>>();
            if coding_segments.is_empty() {
                continue;
            }
            if is_reverse {
                coding_segments.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(b.1.cmp(&a.1)));
            } else {
                coding_segments.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            }

            for (seg_start, seg_end) in coding_segments {
                let seg_len = seg_end.saturating_sub(seg_start);
                if seg_len == 0 {
                    continue;
                }
                let entry_phase = (consumed_cds_bp % 3) as u8;
                let exit_phase = ((consumed_cds_bp + seg_len - 1) % 3) as u8;
                if is_reverse {
                    if seg_end == exon.1 && phases[exon_idx].right_cds_phase.is_none() {
                        phases[exon_idx].right_cds_phase = Some(entry_phase);
                    }
                    if seg_start == exon.0 {
                        phases[exon_idx].left_cds_phase = Some(exit_phase);
                    }
                } else {
                    if seg_start == exon.0 && phases[exon_idx].left_cds_phase.is_none() {
                        phases[exon_idx].left_cds_phase = Some(entry_phase);
                    }
                    if seg_end == exon.1 {
                        phases[exon_idx].right_cds_phase = Some(exit_phase);
                    }
                }
                consumed_cds_bp += seg_len;
            }
        }
        phases
    }

    pub(super) fn build_splicing_expert_view(
        &self,
        seq_id: &str,
        feature_id: usize,
        scope: SplicingScopePreset,
    ) -> Result<SplicingExpertView, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let features = dna.features();
        let target_feature = features.get(feature_id).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "Feature id '{}' was not found in sequence '{}'",
                feature_id, seq_id
            ),
        })?;
        if !Self::is_splicing_seed_feature(target_feature) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Feature '{}' in '{}' is not an mRNA/transcript/ncRNA/misc_RNA/exon/gene/CDS feature and cannot seed a splicing view",
                    feature_id, seq_id
                ),
            });
        }
        let target_group = Self::splicing_group_label(target_feature, feature_id);
        let target_is_reverse = feature_is_reverse(target_feature);
        let restrict_to_target_group = scope.restrict_to_target_group();
        let restrict_to_target_strand = scope.restrict_to_target_strand();

        let mut roi_ranges_0based: Vec<(usize, usize)> = Vec::new();
        for (idx, feature) in features.iter().enumerate() {
            if Self::splicing_group_label(feature, idx) != target_group {
                continue;
            }
            if !Self::is_splicing_transcript_feature(feature) && !Self::is_exon_feature(feature) {
                continue;
            }
            collect_location_ranges_usize(&feature.location, &mut roi_ranges_0based);
        }
        if roi_ranges_0based.is_empty() {
            let (from, to) = target_feature
                .location
                .find_bounds()
                .map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Could not parse target feature range: {e}"),
                })?;
            if from >= 0 && to >= 0 {
                roi_ranges_0based.push((from as usize, to as usize));
            }
        }
        let roi_start_0based = roi_ranges_0based
            .iter()
            .map(|(start, _)| *start)
            .min()
            .unwrap_or(0);
        let roi_end_0based_exclusive = roi_ranges_0based
            .iter()
            .map(|(_, end)| *end)
            .max()
            .unwrap_or(roi_start_0based.saturating_add(1));

        #[derive(Clone)]
        struct TranscriptWork {
            feature_id: usize,
            transcript_id: String,
            label: String,
            is_reverse: bool,
            exon_ranges: Vec<(usize, usize)>,
            exon_cds_phases: Vec<SplicingExonCdsPhase>,
            intron_ranges: Vec<(usize, usize)>,
        }

        let mut transcripts: Vec<TranscriptWork> = Vec::new();
        for (idx, feature) in features.iter().enumerate() {
            if !Self::is_splicing_transcript_feature(feature) {
                continue;
            }
            if restrict_to_target_group && Self::splicing_group_label(feature, idx) != target_group
            {
                continue;
            }
            let is_reverse = feature_is_reverse(feature);
            if restrict_to_target_strand && is_reverse != target_is_reverse {
                continue;
            }
            if !restrict_to_target_group {
                let mut feature_ranges = Vec::<(usize, usize)>::new();
                collect_location_ranges_usize(&feature.location, &mut feature_ranges);
                let overlaps_roi = feature_ranges.iter().any(|(start, end)| {
                    *end > roi_start_0based && *start < roi_end_0based_exclusive
                });
                if !overlaps_roi {
                    continue;
                }
            }
            let mut exon_ranges = vec![];
            collect_location_ranges_usize(&feature.location, &mut exon_ranges);
            if exon_ranges.is_empty() {
                let (from, to) = feature.location.find_bounds().map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Could not parse transcript range: {e}"),
                })?;
                if from >= 0 && to >= 0 {
                    exon_ranges.push((from as usize, to as usize));
                }
            }
            exon_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            exon_ranges.retain(|(start, end)| end > start);
            if exon_ranges.is_empty() {
                continue;
            }
            let mut introns: Vec<(usize, usize)> = vec![];
            for pair in exon_ranges.windows(2) {
                let left = pair[0];
                let right = pair[1];
                if right.0 > left.1 {
                    introns.push((left.1, right.0));
                }
            }
            let exon_cds_phases =
                Self::exon_cds_phases_for_transcript(feature, &exon_ranges, is_reverse);
            transcripts.push(TranscriptWork {
                feature_id: idx,
                transcript_id: Self::feature_transcript_id(feature, idx),
                label: Self::feature_display_label(feature, idx),
                is_reverse,
                exon_ranges,
                exon_cds_phases,
                intron_ranges: introns,
            });
        }

        if transcripts.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No transcript-like RNA features found for splicing group '{}' in '{}'",
                    target_group, seq_id
                ),
            });
        }

        transcripts.sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));

        let mut exon_support: HashMap<(usize, usize), HashSet<usize>> = HashMap::new();
        let mut junction_support: HashMap<(usize, usize), HashSet<usize>> = HashMap::new();
        for transcript in &transcripts {
            for exon in &transcript.exon_ranges {
                exon_support
                    .entry(*exon)
                    .or_default()
                    .insert(transcript.feature_id);
            }
            for intron in &transcript.intron_ranges {
                junction_support
                    .entry(*intron)
                    .or_default()
                    .insert(transcript.feature_id);
            }
        }

        let mut unique_exons: Vec<(usize, usize)> = exon_support.keys().copied().collect();
        unique_exons.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let region_start = unique_exons
            .iter()
            .map(|(start, _)| *start)
            .min()
            .unwrap_or(0);
        let region_end = unique_exons
            .iter()
            .map(|(_, end)| *end)
            .max()
            .unwrap_or(region_start + 1);

        let transcript_count = transcripts.len();
        let unique_exon_summaries: Vec<SplicingExonSummary> = unique_exons
            .iter()
            .map(|(start, end)| {
                let support = exon_support
                    .get(&(*start, *end))
                    .map(|set| set.len())
                    .unwrap_or(0);
                SplicingExonSummary {
                    start_1based: start + 1,
                    end_1based: *end,
                    support_transcript_count: support,
                    constitutive: support == transcript_count,
                }
            })
            .collect();

        let transcript_lanes: Vec<SplicingTranscriptLane> = transcripts
            .iter()
            .map(|transcript| SplicingTranscriptLane {
                transcript_feature_id: transcript.feature_id,
                transcript_id: transcript.transcript_id.clone(),
                label: transcript.label.clone(),
                strand: if transcript.is_reverse {
                    "-".to_string()
                } else {
                    "+".to_string()
                },
                exons: Self::range_vec_to_splicing(transcript.exon_ranges.clone()),
                exon_cds_phases: transcript.exon_cds_phases.clone(),
                introns: Self::range_vec_to_splicing(transcript.intron_ranges.clone()),
                has_target_feature: transcript.feature_id == feature_id,
            })
            .collect();

        let matrix_rows: Vec<SplicingMatrixRow> = transcripts
            .iter()
            .map(|transcript| {
                let set = transcript
                    .exon_ranges
                    .iter()
                    .copied()
                    .collect::<HashSet<(usize, usize)>>();
                SplicingMatrixRow {
                    transcript_feature_id: transcript.feature_id,
                    transcript_id: transcript.transcript_id.clone(),
                    label: transcript.label.clone(),
                    exon_presence: unique_exons.iter().map(|exon| set.contains(exon)).collect(),
                }
            })
            .collect();

        let mut boundaries = Vec::new();
        for transcript in &transcripts {
            boundaries.extend(Self::splice_boundary_markers_for_introns(
                dna,
                transcript.feature_id,
                &transcript.transcript_id,
                transcript.is_reverse,
                &transcript.intron_ranges,
            ));
        }
        boundaries.sort_by(|a, b| {
            a.position_1based
                .cmp(&b.position_1based)
                .then_with(|| a.side.cmp(&b.side))
                .then_with(|| a.transcript_feature_id.cmp(&b.transcript_feature_id))
        });

        let mut junctions: Vec<SplicingJunctionArc> = junction_support
            .into_iter()
            .map(|((donor0, acceptor0), ids)| {
                let mut transcript_feature_ids = ids.into_iter().collect::<Vec<_>>();
                transcript_feature_ids.sort_unstable();
                SplicingJunctionArc {
                    donor_1based: donor0 + 1,
                    acceptor_1based: acceptor0,
                    support_transcript_count: transcript_feature_ids.len(),
                    transcript_feature_ids,
                }
            })
            .collect();
        junctions.sort_by(|a, b| {
            a.donor_1based
                .cmp(&b.donor_1based)
                .then(a.acceptor_1based.cmp(&b.acceptor_1based))
        });

        let mut exon_skip_details = Vec::new();
        let mut intron_retention_details = Vec::new();
        let mut alt5_details = Vec::new();
        let mut alt3_details = Vec::new();
        let mut mex_details = Vec::new();

        let exon_set_by_transcript: HashMap<usize, HashSet<(usize, usize)>> = transcripts
            .iter()
            .map(|transcript| {
                (
                    transcript.feature_id,
                    transcript.exon_ranges.iter().copied().collect(),
                )
            })
            .collect();

        let mut alt5_groups = 0usize;
        let mut alt3_groups = 0usize;
        let mut grouped_by_start: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
        let mut grouped_by_end: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
        for exon in &unique_exons {
            grouped_by_start.entry(exon.0).or_default().push(*exon);
            grouped_by_end.entry(exon.1).or_default().push(*exon);
        }
        if target_is_reverse {
            for (start, group) in grouped_by_start {
                let distinct = group
                    .iter()
                    .map(|(_, end)| *end)
                    .collect::<HashSet<_>>()
                    .len();
                if distinct > 1 {
                    alt5_groups += 1;
                    alt5_details.push(format!(
                        "shared start={} with {} alternative 5' exon end(s)",
                        start + 1,
                        distinct
                    ));
                }
            }
            for (end, group) in grouped_by_end {
                let distinct = group
                    .iter()
                    .map(|(start, _)| *start)
                    .collect::<HashSet<_>>()
                    .len();
                if distinct > 1 {
                    alt3_groups += 1;
                    alt3_details.push(format!(
                        "shared end={} with {} alternative 3' exon start(s)",
                        end, distinct
                    ));
                }
            }
        } else {
            for (start, group) in grouped_by_start {
                let distinct = group
                    .iter()
                    .map(|(_, end)| *end)
                    .collect::<HashSet<_>>()
                    .len();
                if distinct > 1 {
                    alt3_groups += 1;
                    alt3_details.push(format!(
                        "shared start={} with {} alternative 3' exon end(s)",
                        start + 1,
                        distinct
                    ));
                }
            }
            for (end, group) in grouped_by_end {
                let distinct = group
                    .iter()
                    .map(|(start, _)| *start)
                    .collect::<HashSet<_>>()
                    .len();
                if distinct > 1 {
                    alt5_groups += 1;
                    alt5_details.push(format!(
                        "shared end={} with {} alternative 5' exon start(s)",
                        end, distinct
                    ));
                }
            }
        }

        let mut exon_skip_count = 0usize;
        for transcript in &transcripts {
            for (intron_start, intron_end) in &transcript.intron_ranges {
                let skipped = unique_exons
                    .iter()
                    .filter(|(exon_start, exon_end)| {
                        *exon_start >= *intron_start
                            && *exon_end <= *intron_end
                            && !transcript
                                .exon_ranges
                                .iter()
                                .any(|own| own.0 == *exon_start && own.1 == *exon_end)
                    })
                    .count();
                if skipped > 0 {
                    exon_skip_count += skipped;
                    exon_skip_details.push(format!(
                        "{} junction {}..{} skips {} exon(s)",
                        transcript.transcript_id,
                        intron_start + 1,
                        intron_end,
                        skipped
                    ));
                }
            }
        }

        let mut intron_retention_count = 0usize;
        for transcript in &transcripts {
            for (intron_start, intron_end) in &transcript.intron_ranges {
                let retained_by = transcripts.iter().find_map(|other| {
                    if other.feature_id == transcript.feature_id {
                        return None;
                    }
                    if other.exon_ranges.iter().any(|(exon_start, exon_end)| {
                        exon_start <= intron_start && exon_end >= intron_end
                    }) {
                        Some(other.transcript_id.clone())
                    } else {
                        None
                    }
                });
                if let Some(other_id) = retained_by {
                    intron_retention_count += 1;
                    intron_retention_details.push(format!(
                        "{} intron {}..{} retained in {}",
                        transcript.transcript_id,
                        intron_start + 1,
                        intron_end,
                        other_id
                    ));
                }
            }
        }

        let mut mex_count = 0usize;
        for left_idx in 0..unique_exons.len() {
            for right_idx in (left_idx + 1)..unique_exons.len() {
                let left = unique_exons[left_idx];
                let right = unique_exons[right_idx];
                if left.1 <= right.0 || right.1 <= left.0 {
                    continue;
                }
                let left_support = exon_support.get(&left).cloned().unwrap_or_default();
                let right_support = exon_support.get(&right).cloned().unwrap_or_default();
                if left_support.is_empty()
                    || right_support.is_empty()
                    || !left_support.is_disjoint(&right_support)
                {
                    continue;
                }
                let co_occurs = transcripts.iter().any(|transcript| {
                    let set = exon_set_by_transcript
                        .get(&transcript.feature_id)
                        .cloned()
                        .unwrap_or_default();
                    set.contains(&left) && set.contains(&right)
                });
                if !co_occurs {
                    mex_count += 1;
                    mex_details.push(format!(
                        "{}..{} and {}..{} appear mutually exclusive",
                        left.0 + 1,
                        left.1,
                        right.0 + 1,
                        right.1
                    ));
                }
            }
        }

        let alt_exon_count = unique_exon_summaries
            .iter()
            .filter(|exon| !exon.constitutive)
            .count();

        let events = vec![
            SplicingEventSummary {
                event_type: "alternative_exon".to_string(),
                count: alt_exon_count,
                details: unique_exon_summaries
                    .iter()
                    .filter(|exon| !exon.constitutive)
                    .take(6)
                    .map(|exon| {
                        format!(
                            "{}..{} support={}/{}",
                            exon.start_1based,
                            exon.end_1based,
                            exon.support_transcript_count,
                            transcript_count
                        )
                    })
                    .collect(),
            },
            SplicingEventSummary {
                event_type: "alternative_5prime".to_string(),
                count: alt5_groups,
                details: alt5_details.into_iter().take(6).collect(),
            },
            SplicingEventSummary {
                event_type: "alternative_3prime".to_string(),
                count: alt3_groups,
                details: alt3_details.into_iter().take(6).collect(),
            },
            SplicingEventSummary {
                event_type: "exon_skipping".to_string(),
                count: exon_skip_count,
                details: exon_skip_details.into_iter().take(6).collect(),
            },
            SplicingEventSummary {
                event_type: "intron_retention".to_string(),
                count: intron_retention_count,
                details: intron_retention_details.into_iter().take(6).collect(),
            },
            SplicingEventSummary {
                event_type: "mutually_exclusive_exons".to_string(),
                count: mex_count,
                details: mex_details.into_iter().take(6).collect(),
            },
        ];

        Ok(SplicingExpertView {
            seq_id: seq_id.to_string(),
            target_feature_id: feature_id,
            group_label: target_group,
            strand: if target_is_reverse {
                "-".to_string()
            } else {
                "+".to_string()
            },
            region_start_1based: region_start + 1,
            region_end_1based: region_end,
            transcript_count,
            unique_exon_count: unique_exons.len(),
            instruction: SPLICING_EXPERT_INSTRUCTION.to_string(),
            transcripts: transcript_lanes,
            unique_exons: unique_exon_summaries,
            matrix_rows,
            boundaries,
            junctions,
            events,
        })
    }

    pub fn inspect_feature_expert(
        &self,
        seq_id: &str,
        target: &FeatureExpertTarget,
    ) -> Result<FeatureExpertView, EngineError> {
        match target {
            FeatureExpertTarget::TfbsFeature { feature_id } => self
                .build_tfbs_expert_view(seq_id, *feature_id)
                .map(FeatureExpertView::Tfbs),
            FeatureExpertTarget::RestrictionSite {
                cut_pos_1based,
                enzyme,
                recognition_start_1based,
                recognition_end_1based,
            } => self
                .build_restriction_site_expert_view(
                    seq_id,
                    *cut_pos_1based,
                    enzyme.as_deref(),
                    *recognition_start_1based,
                    *recognition_end_1based,
                )
                .map(FeatureExpertView::RestrictionSite),
            FeatureExpertTarget::SplicingFeature { feature_id, scope } => self
                .build_splicing_expert_view(seq_id, *feature_id, *scope)
                .map(FeatureExpertView::Splicing),
            FeatureExpertTarget::IsoformArchitecture { panel_id } => self
                .build_isoform_architecture_expert_view(seq_id, panel_id)
                .map(FeatureExpertView::IsoformArchitecture),
        }
    }

    pub fn render_feature_expert_svg_to_path(
        &self,
        seq_id: &str,
        target: &FeatureExpertTarget,
        path: &str,
    ) -> Result<FeatureExpertView, EngineError> {
        let view = self.inspect_feature_expert(seq_id, target)?;
        let svg = render_feature_expert_svg(&view);
        std::fs::write(path, svg).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write feature-expert SVG to '{path}': {e}"),
        })?;
        Ok(view)
    }

    pub fn render_isoform_architecture_svg_to_path(
        &self,
        seq_id: &str,
        panel_id: &str,
        path: &str,
    ) -> Result<IsoformArchitectureExpertView, EngineError> {
        let target = FeatureExpertTarget::IsoformArchitecture {
            panel_id: panel_id.to_string(),
        };
        let view = self.render_feature_expert_svg_to_path(seq_id, &target, path)?;
        match view {
            FeatureExpertView::IsoformArchitecture(isoform) => Ok(isoform),
            _ => Err(EngineError {
                code: ErrorCode::Internal,
                message: "Unexpected expert-view payload while rendering isoform architecture SVG"
                    .to_string(),
            }),
        }
    }
}
