//! Candidate-set summary, metric computation, and expression-evaluation helpers.
//!
//! Keep reusable candidate inspection/reporting code here so candidate-related
//! operations and adapters share one scoring vocabulary.
//!
//! Look here for:
//! - candidate-set paging and summary views
//! - metric-name discovery and weighted-objective helpers
//! - reusable sort/filter/scoring logic that should stay independent from
//!   shell/GUI command routing

use super::*;

impl GentleEngine {
    pub(super) fn metric_names_for_candidate_set(set: &CandidateSet) -> Vec<String> {
        let mut names = BTreeSet::new();
        for candidate in &set.candidates {
            for metric in candidate.metrics.keys() {
                names.insert(metric.to_string());
            }
        }
        names.into_iter().collect()
    }

    pub(super) fn candidate_key(record: &CandidateRecord) -> String {
        format!(
            "{}:{}:{}:{}",
            record.seq_id, record.start_0based, record.end_0based, record.sequence
        )
    }

    /// Summarize all stored candidate sets in deterministic name order.
    pub fn list_candidate_sets(&self) -> Vec<CandidateSetSummary> {
        let store = self.read_candidate_store();
        let mut set_names: Vec<String> = store.sets.keys().cloned().collect();
        set_names.sort();
        set_names
            .iter()
            .filter_map(|name| store.sets.get(name))
            .map(|set| CandidateSetSummary {
                name: set.name.clone(),
                created_at_unix_ms: set.created_at_unix_ms,
                source_seq_ids: set.source_seq_ids.clone(),
                candidate_count: set.candidates.len(),
                metrics: Self::metric_names_for_candidate_set(set),
            })
            .collect()
    }

    /// Return one paged candidate-set slice plus total-count pagination data.
    pub fn inspect_candidate_set_page(
        &self,
        set_name: &str,
        limit: usize,
        offset: usize,
    ) -> Result<(CandidateSet, usize, usize), EngineError> {
        if limit == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Candidate page limit must be >= 1".to_string(),
            });
        }
        let set_name = Self::normalize_candidate_set_name(set_name)?;
        let store = self.read_candidate_store();
        let mut set = store
            .sets
            .get(&set_name)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Candidate set '{}' not found", set_name),
            })?;
        let total = set.candidates.len();
        let clamped_offset = offset.min(total);
        set.candidates = set
            .candidates
            .into_iter()
            .skip(clamped_offset)
            .take(limit)
            .collect();
        Ok((set, total, clamped_offset))
    }

    /// Summarize which metric keys are present across one candidate set.
    pub fn list_candidate_set_metrics(
        &self,
        set_name: &str,
    ) -> Result<Vec<CandidateMetricSummary>, EngineError> {
        let set_name = Self::normalize_candidate_set_name(set_name)?;
        let store = self.read_candidate_store();
        let set = store.sets.get(&set_name).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!("Candidate set '{}' not found", set_name),
        })?;
        let mut counts: HashMap<String, usize> = HashMap::new();
        for candidate in &set.candidates {
            for metric in candidate.metrics.keys() {
                *counts.entry(metric.to_string()).or_insert(0) += 1;
            }
        }
        let mut summaries = counts
            .into_iter()
            .map(|(metric, present_in_candidates)| CandidateMetricSummary {
                metric,
                present_in_candidates,
                missing_in_candidates: set.candidates.len().saturating_sub(present_in_candidates),
            })
            .collect::<Vec<_>>();
        summaries.sort_by(|a, b| a.metric.cmp(&b.metric));
        Ok(summaries)
    }

    pub(super) fn sequence_base_counts(sequence: &[u8]) -> (usize, usize, usize, usize, usize) {
        let mut a = 0usize;
        let mut c = 0usize;
        let mut g = 0usize;
        let mut t = 0usize;
        for b in sequence {
            match b.to_ascii_uppercase() {
                b'A' => a += 1,
                b'C' => c += 1,
                b'G' => g += 1,
                b'T' => t += 1,
                _ => {}
            }
        }
        let canonical = a + c + g + t;
        (a, c, g, t, canonical)
    }

    pub(super) fn approximate_molecular_weight_da(sequence: &[u8]) -> f64 {
        let mut total = 0.0;
        for b in sequence {
            total += match b.to_ascii_uppercase() {
                b'A' => 313.21,
                b'C' => 289.18,
                b'G' => 329.21,
                b'T' => 304.20,
                _ => 0.0,
            };
        }
        total
    }

    pub(super) fn feature_labels_upper(feature: &gb_io::seq::Feature) -> Vec<String> {
        let mut labels = vec![];
        for key in [
            "label",
            "gene",
            "gene_name",
            "locus_tag",
            "product",
            "standard_name",
            "note",
        ] {
            for value in feature.qualifier_values(key) {
                let trimmed = value.trim();
                if !trimmed.is_empty() {
                    labels.push(trimmed.to_ascii_uppercase());
                }
            }
        }
        labels
    }

    pub(super) fn collect_location_segments(
        location: &gb_io::seq::Location,
        seq_len: usize,
        reverse: bool,
        out: &mut Vec<FeatureLocationSegment>,
    ) {
        use gb_io::seq::Location;
        match location {
            Location::Range((a, _), (b, _)) => {
                let raw_start = (*a).min(*b);
                let raw_end = (*a).max(*b);
                let raw_end_exclusive = if raw_end == raw_start {
                    raw_end.saturating_add(1)
                } else {
                    raw_end
                };
                if raw_end_exclusive <= 0 || raw_start >= seq_len as i64 {
                    return;
                }
                let start = raw_start.max(0).min(seq_len as i64) as usize;
                let end_0based = raw_end_exclusive.max(0).min(seq_len as i64) as usize;
                if end_0based <= start {
                    return;
                }
                out.push(FeatureLocationSegment {
                    start_0based: start,
                    end_0based,
                    strand: Some(if reverse { '-' } else { '+' }),
                });
            }
            Location::Between(a, b) => {
                let raw_start = (*a).min(*b);
                let raw_end_exclusive = (*a).max(*b).saturating_add(1);
                if raw_end_exclusive <= 0 || raw_start >= seq_len as i64 {
                    return;
                }
                let start = raw_start.max(0).min(seq_len as i64) as usize;
                let end_0based = raw_end_exclusive.max(0).min(seq_len as i64) as usize;
                if end_0based <= start {
                    return;
                }
                out.push(FeatureLocationSegment {
                    start_0based: start,
                    end_0based,
                    strand: Some(if reverse { '-' } else { '+' }),
                });
            }
            Location::Complement(inner) => {
                Self::collect_location_segments(inner, seq_len, !reverse, out);
            }
            Location::Join(parts)
            | Location::Order(parts)
            | Location::Bond(parts)
            | Location::OneOf(parts) => {
                for part in parts {
                    Self::collect_location_segments(part, seq_len, reverse, out);
                }
            }
            Location::External(_, maybe_inner) => {
                if let Some(inner) = maybe_inner.as_deref() {
                    Self::collect_location_segments(inner, seq_len, reverse, out);
                }
            }
            Location::Gap(_) => {}
        }
    }

    pub(super) fn feature_segment_boundary_positions(
        segment: &FeatureLocationSegment,
        boundary_mode: CandidateFeatureBoundaryMode,
    ) -> Vec<usize> {
        if segment.end_0based <= segment.start_0based {
            return vec![];
        }
        let start = segment.start_0based;
        let end = segment.end_0based.saturating_sub(1);
        let mut out = vec![];
        match boundary_mode {
            CandidateFeatureBoundaryMode::Any => {
                out.push(start);
                out.push(end);
            }
            CandidateFeatureBoundaryMode::Start => out.push(start),
            CandidateFeatureBoundaryMode::End => out.push(end),
            CandidateFeatureBoundaryMode::FivePrime => match segment.strand {
                Some('+') => out.push(start),
                Some('-') => out.push(end),
                _ => {
                    out.push(start);
                    out.push(end);
                }
            },
            CandidateFeatureBoundaryMode::ThreePrime => match segment.strand {
                Some('+') => out.push(end),
                Some('-') => out.push(start),
                _ => {
                    out.push(start);
                    out.push(end);
                }
            },
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    pub(super) fn collect_feature_distance_targets(
        dna: &DNAsequence,
        geometry_mode: CandidateFeatureGeometryMode,
        boundary_mode: CandidateFeatureBoundaryMode,
    ) -> Vec<FeatureDistanceTarget> {
        let mut out = vec![];
        for (feature_index, feature) in dna.features().iter().enumerate() {
            if feature.kind.to_string().eq_ignore_ascii_case("SOURCE") {
                continue;
            }
            let kind_upper = feature.kind.to_string().to_ascii_uppercase();
            let labels_upper = Self::feature_labels_upper(feature);
            let mut segments = vec![];
            Self::collect_location_segments(&feature.location, dna.len(), false, &mut segments);
            if segments.is_empty() {
                let Ok((from, to)) = feature.location.find_bounds() else {
                    continue;
                };
                if from < 0 || to <= 0 {
                    continue;
                }
                let start = from.min(to).max(0).min(dna.len() as i64) as usize;
                let raw_end = from.max(to);
                let raw_end_exclusive = if raw_end == from.min(to) {
                    raw_end.saturating_add(1)
                } else {
                    raw_end
                };
                let end_0based = raw_end_exclusive.max(0).min(dna.len() as i64) as usize;
                if end_0based <= start {
                    continue;
                }
                segments.push(FeatureLocationSegment {
                    start_0based: start,
                    end_0based,
                    strand: None,
                });
            }
            if segments.is_empty() {
                continue;
            }
            let feature_strand = {
                let mut strand = None;
                let mut ambiguous = false;
                for segment in &segments {
                    match (strand, segment.strand) {
                        (_, None) => {
                            ambiguous = true;
                            break;
                        }
                        (None, Some(s)) => strand = Some(s),
                        (Some(prev), Some(s)) if prev == s => {}
                        (Some(_), Some(_)) => {
                            ambiguous = true;
                            break;
                        }
                    }
                }
                if ambiguous { None } else { strand }
            };
            match geometry_mode {
                CandidateFeatureGeometryMode::FeatureSpan => {
                    let start = segments
                        .iter()
                        .map(|segment| segment.start_0based)
                        .min()
                        .unwrap_or(0);
                    let end_0based = segments
                        .iter()
                        .map(|segment| segment.end_0based)
                        .max()
                        .unwrap_or(start);
                    if end_0based <= start {
                        continue;
                    }
                    out.push(FeatureDistanceTarget {
                        feature_index,
                        kind_upper: kind_upper.clone(),
                        labels_upper: labels_upper.clone(),
                        start_0based: start,
                        end_0based,
                        strand: feature_strand,
                    });
                }
                CandidateFeatureGeometryMode::FeatureParts => {
                    for segment in segments {
                        if segment.end_0based <= segment.start_0based {
                            continue;
                        }
                        out.push(FeatureDistanceTarget {
                            feature_index,
                            kind_upper: kind_upper.clone(),
                            labels_upper: labels_upper.clone(),
                            start_0based: segment.start_0based,
                            end_0based: segment.end_0based,
                            strand: segment.strand,
                        });
                    }
                }
                CandidateFeatureGeometryMode::FeatureBoundaries => {
                    for segment in &segments {
                        for point in
                            Self::feature_segment_boundary_positions(segment, boundary_mode)
                        {
                            if point >= dna.len() {
                                continue;
                            }
                            out.push(FeatureDistanceTarget {
                                feature_index,
                                kind_upper: kind_upper.clone(),
                                labels_upper: labels_upper.clone(),
                                start_0based: point,
                                end_0based: point.saturating_add(1).min(dna.len()),
                                strand: segment.strand,
                            });
                        }
                    }
                }
            }
        }
        out
    }

    pub(super) fn compile_optional_regex(
        pattern: &Option<String>,
        option_name: &str,
    ) -> Result<Option<Regex>, EngineError> {
        let Some(raw) = pattern else {
            return Ok(None);
        };
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Ok(None);
        }
        RegexBuilder::new(trimmed)
            .case_insensitive(true)
            .build()
            .map(Some)
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Invalid {option_name} regex '{}': {}", trimmed, e),
            })
    }

    pub(super) fn feature_matches_filter(
        feature: &FeatureDistanceTarget,
        kind_filter_upper: &[String],
        label_regex: Option<&Regex>,
    ) -> bool {
        let kind_ok = kind_filter_upper.is_empty()
            || kind_filter_upper
                .iter()
                .any(|kind| feature.kind_upper == *kind);
        if !kind_ok {
            return false;
        }
        if let Some(regex) = label_regex {
            feature
                .labels_upper
                .iter()
                .any(|label| regex.is_match(label))
        } else {
            true
        }
    }

    pub(super) fn feature_matches_strand_relation(
        feature: &FeatureDistanceTarget,
        strand_relation: CandidateFeatureStrandRelation,
    ) -> bool {
        match strand_relation {
            CandidateFeatureStrandRelation::Any => true,
            CandidateFeatureStrandRelation::Same => feature.strand == Some('+'),
            CandidateFeatureStrandRelation::Opposite => feature.strand == Some('-'),
        }
    }

    pub(super) fn interval_distance(
        left_start: usize,
        left_end: usize,
        right_start: usize,
        right_end: usize,
    ) -> usize {
        if left_end <= right_start {
            right_start.saturating_sub(left_end)
        } else if right_end <= left_start {
            left_start.saturating_sub(right_end)
        } else {
            0
        }
    }

    pub(super) fn boundary_distance(
        candidate_start: usize,
        candidate_end: usize,
        boundary_pos: usize,
    ) -> usize {
        if boundary_pos < candidate_start {
            candidate_start.saturating_sub(boundary_pos)
        } else if boundary_pos > candidate_end {
            boundary_pos.saturating_sub(candidate_end)
        } else {
            0
        }
    }

    pub(super) fn nearest_feature_distance(
        candidate_start: usize,
        candidate_end: usize,
        features: &[FeatureDistanceTarget],
        kind_filter_upper: &[String],
        label_regex: Option<&Regex>,
        strand_relation: CandidateFeatureStrandRelation,
    ) -> Option<usize> {
        features
            .iter()
            .filter(|feature| Self::feature_matches_filter(feature, kind_filter_upper, label_regex))
            .filter(|feature| Self::feature_matches_strand_relation(feature, strand_relation))
            .map(|feature| {
                Self::interval_distance(
                    candidate_start,
                    candidate_end,
                    feature.start_0based,
                    feature.end_0based,
                )
            })
            .min()
    }

    pub(super) fn matching_feature_count(
        features: &[FeatureDistanceTarget],
        kind_filter_upper: &[String],
        label_regex: Option<&Regex>,
        strand_relation: CandidateFeatureStrandRelation,
    ) -> usize {
        features
            .iter()
            .filter(|feature| Self::feature_matches_filter(feature, kind_filter_upper, label_regex))
            .filter(|feature| Self::feature_matches_strand_relation(feature, strand_relation))
            .map(|feature| feature.feature_index)
            .collect::<HashSet<_>>()
            .len()
    }

    pub(super) fn compute_candidate_metrics(
        sequence: &[u8],
        start_0based: usize,
        end_0based: usize,
        source_len: usize,
        nearest_feature_distance_bp: Option<usize>,
    ) -> HashMap<String, f64> {
        let mut metrics = HashMap::new();
        let length_bp = end_0based.saturating_sub(start_0based);
        let (count_a, count_c, count_g, count_t, canonical) = Self::sequence_base_counts(sequence);
        metrics.insert("length_bp".to_string(), length_bp as f64);
        metrics.insert(
            "molecular_weight_da".to_string(),
            Self::approximate_molecular_weight_da(sequence),
        );
        metrics.insert("distance_to_seq_start_bp".to_string(), start_0based as f64);
        metrics.insert(
            "distance_to_seq_end_bp".to_string(),
            source_len.saturating_sub(end_0based) as f64,
        );
        metrics.insert("count_a".to_string(), count_a as f64);
        metrics.insert("count_c".to_string(), count_c as f64);
        metrics.insert("count_g".to_string(), count_g as f64);
        metrics.insert("count_t".to_string(), count_t as f64);
        if canonical > 0 {
            let canonical_f = canonical as f64;
            metrics.insert(
                "gc_fraction".to_string(),
                (count_g + count_c) as f64 / canonical_f,
            );
            metrics.insert(
                "at_fraction".to_string(),
                (count_a + count_t) as f64 / canonical_f,
            );
            metrics.insert("a_fraction".to_string(), count_a as f64 / canonical_f);
            metrics.insert("c_fraction".to_string(), count_c as f64 / canonical_f);
            metrics.insert("g_fraction".to_string(), count_g as f64 / canonical_f);
            metrics.insert("t_fraction".to_string(), count_t as f64 / canonical_f);
        }
        if count_t > 0 {
            metrics.insert("at_ratio".to_string(), count_a as f64 / count_t as f64);
        }
        if count_c > 0 {
            metrics.insert("gc_ratio".to_string(), count_g as f64 / count_c as f64);
        }
        if let Some(distance) = nearest_feature_distance_bp {
            metrics.insert(
                "distance_to_nearest_feature_bp".to_string(),
                distance as f64,
            );
        }
        metrics.retain(|_, value| value.is_finite());
        metrics
    }

    pub(super) fn quantile_threshold(values: &[f64], quantile: f64) -> Option<f64> {
        if values.is_empty() || !quantile.is_finite() || !(0.0..=1.0).contains(&quantile) {
            return None;
        }
        let mut sorted = values.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        let max_idx = sorted.len().saturating_sub(1);
        let idx = ((max_idx as f64) * quantile).round() as usize;
        sorted.get(idx.min(max_idx)).copied()
    }

    pub(super) fn tokenize_expression(
        expression: &str,
    ) -> Result<Vec<ExpressionToken>, EngineError> {
        let bytes = expression.as_bytes();
        let mut idx = 0usize;
        let mut tokens = vec![];
        while idx < bytes.len() {
            let b = bytes[idx];
            if b.is_ascii_whitespace() {
                idx += 1;
                continue;
            }
            match b {
                b'+' => {
                    tokens.push(ExpressionToken::Plus);
                    idx += 1;
                }
                b'-' => {
                    tokens.push(ExpressionToken::Minus);
                    idx += 1;
                }
                b'*' => {
                    tokens.push(ExpressionToken::Star);
                    idx += 1;
                }
                b'/' => {
                    tokens.push(ExpressionToken::Slash);
                    idx += 1;
                }
                b'(' => {
                    tokens.push(ExpressionToken::LParen);
                    idx += 1;
                }
                b')' => {
                    tokens.push(ExpressionToken::RParen);
                    idx += 1;
                }
                b',' => {
                    tokens.push(ExpressionToken::Comma);
                    idx += 1;
                }
                _ if b.is_ascii_digit() || b == b'.' => {
                    let start = idx;
                    idx += 1;
                    while idx < bytes.len()
                        && (bytes[idx].is_ascii_digit()
                            || bytes[idx] == b'.'
                            || matches!(bytes[idx], b'e' | b'E' | b'+' | b'-'))
                    {
                        if matches!(bytes[idx], b'+' | b'-') {
                            let prev = bytes[idx.saturating_sub(1)];
                            if !matches!(prev, b'e' | b'E') {
                                break;
                            }
                        }
                        idx += 1;
                    }
                    let raw = &expression[start..idx];
                    let value = raw.parse::<f64>().map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Invalid numeric literal '{}': {}", raw, e),
                    })?;
                    tokens.push(ExpressionToken::Number(value));
                }
                _ if b.is_ascii_alphabetic() || b == b'_' => {
                    let start = idx;
                    idx += 1;
                    while idx < bytes.len()
                        && (bytes[idx].is_ascii_alphanumeric() || bytes[idx] == b'_')
                    {
                        idx += 1;
                    }
                    tokens.push(ExpressionToken::Ident(expression[start..idx].to_string()));
                }
                _ => {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Unsupported character '{}' in expression", b as char),
                    });
                }
            }
        }
        if tokens.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Expression is empty".to_string(),
            });
        }
        Ok(tokens)
    }

    pub(super) fn parse_metric_expression(expression: &str) -> Result<MetricExpr, EngineError> {
        let tokens = Self::tokenize_expression(expression)?;
        MetricExpressionParser::new(tokens)
            .parse()
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Invalid score expression: {e}"),
            })
    }

    pub(super) fn evaluate_metric_expression(
        expr: &MetricExpr,
        metrics: &HashMap<String, f64>,
    ) -> Result<f64, EngineError> {
        match expr {
            MetricExpr::Number(value) => Ok(*value),
            MetricExpr::Variable(name) => {
                if let Some(value) = metrics.get(name) {
                    return Ok(*value);
                }
                let normalized = Self::normalize_metric_name(name);
                metrics
                    .get(&normalized)
                    .copied()
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Expression references unknown metric '{}'", name),
                    })
            }
            MetricExpr::UnaryMinus(inner) => Ok(-Self::evaluate_metric_expression(inner, metrics)?),
            MetricExpr::Binary { op, left, right } => {
                let lhs = Self::evaluate_metric_expression(left, metrics)?;
                let rhs = Self::evaluate_metric_expression(right, metrics)?;
                let value = match op {
                    ExpressionBinaryOp::Add => lhs + rhs,
                    ExpressionBinaryOp::Subtract => lhs - rhs,
                    ExpressionBinaryOp::Multiply => lhs * rhs,
                    ExpressionBinaryOp::Divide => {
                        if rhs.abs() < f64::EPSILON {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Expression division by zero".to_string(),
                            });
                        }
                        lhs / rhs
                    }
                };
                if value.is_finite() {
                    Ok(value)
                } else {
                    Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Expression produced a non-finite value".to_string(),
                    })
                }
            }
            MetricExpr::Function { name, args } => {
                let normalized = name.trim().to_ascii_lowercase();
                let value = match normalized.as_str() {
                    "abs" => {
                        if args.len() != 1 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function abs() expects exactly 1 argument".to_string(),
                            });
                        }
                        Self::evaluate_metric_expression(&args[0], metrics)?.abs()
                    }
                    "sqrt" => {
                        if args.len() != 1 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function sqrt() expects exactly 1 argument".to_string(),
                            });
                        }
                        let x = Self::evaluate_metric_expression(&args[0], metrics)?;
                        if x < 0.0 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function sqrt() requires non-negative input".to_string(),
                            });
                        }
                        x.sqrt()
                    }
                    "log" => {
                        if args.len() != 1 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function log() expects exactly 1 argument".to_string(),
                            });
                        }
                        let x = Self::evaluate_metric_expression(&args[0], metrics)?;
                        if x <= 0.0 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function log() requires positive input".to_string(),
                            });
                        }
                        x.ln()
                    }
                    "exp" => {
                        if args.len() != 1 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function exp() expects exactly 1 argument".to_string(),
                            });
                        }
                        Self::evaluate_metric_expression(&args[0], metrics)?.exp()
                    }
                    "min" => {
                        if args.len() != 2 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function min() expects exactly 2 arguments".to_string(),
                            });
                        }
                        let a = Self::evaluate_metric_expression(&args[0], metrics)?;
                        let b = Self::evaluate_metric_expression(&args[1], metrics)?;
                        a.min(b)
                    }
                    "max" => {
                        if args.len() != 2 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function max() expects exactly 2 arguments".to_string(),
                            });
                        }
                        let a = Self::evaluate_metric_expression(&args[0], metrics)?;
                        let b = Self::evaluate_metric_expression(&args[1], metrics)?;
                        a.max(b)
                    }
                    "pow" => {
                        if args.len() != 2 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function pow() expects exactly 2 arguments".to_string(),
                            });
                        }
                        let a = Self::evaluate_metric_expression(&args[0], metrics)?;
                        let b = Self::evaluate_metric_expression(&args[1], metrics)?;
                        a.powf(b)
                    }
                    "clamp" => {
                        if args.len() != 3 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function clamp() expects exactly 3 arguments".to_string(),
                            });
                        }
                        let x = Self::evaluate_metric_expression(&args[0], metrics)?;
                        let lo = Self::evaluate_metric_expression(&args[1], metrics)?;
                        let hi = Self::evaluate_metric_expression(&args[2], metrics)?;
                        if lo > hi {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function clamp() requires lo <= hi".to_string(),
                            });
                        }
                        x.max(lo).min(hi)
                    }
                    _ => {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("Unknown expression function '{}'", name),
                        });
                    }
                };
                if value.is_finite() {
                    Ok(value)
                } else {
                    Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Function '{}' produced a non-finite value", name),
                    })
                }
            }
        }
    }
}
