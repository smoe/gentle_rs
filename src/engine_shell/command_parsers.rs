//! Extracted command parser helpers for candidate/guide/macro/routine/planning
//! families.
//!
//! Keep large command-family token parsing here so `engine_shell.rs` can remain
//! the readable public contract surface while still routing every family
//! through one shared parser.
//!
//! Look here for:
//! - token-to-`ShellCommand` parsing for the larger command families
//! - family-specific option validation that should stay shared across GUI shell
//!   and CLI shell mode
//! - places to extend when a shell family gets too large for `engine_shell.rs`

use super::*;
use crate::engine::{
    TfbsScoreTrackCorrelationMetric, TfbsScoreTrackCorrelationSignalSource,
    TfbsScoreTrackValueKind, TfbsTrackSimilarityRankingMetric,
};

pub(super) fn parse_containers_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("containers requires a subcommand: set-exclusive".to_string());
    }
    match tokens[1].as_str() {
        "set-exclusive" => {
            if tokens.len() != 4 {
                return Err("containers set-exclusive requires CONTAINER_ID true|false".to_string());
            }
            let container_id = tokens[2].trim().to_string();
            if container_id.is_empty() {
                return Err(
                    "containers set-exclusive requires a non-empty CONTAINER_ID".to_string()
                );
            }
            let exclusive = match tokens[3].trim().to_ascii_lowercase().as_str() {
                "true" | "yes" | "on" | "exclusive" | "declared_only" | "declared-only" => true,
                "false" | "no" | "off" | "subset" | "known_subset" | "known-subset" => false,
                other => {
                    return Err(format!(
                        "Unsupported exclusivity value '{other}', expected true|false"
                    ));
                }
            };
            Ok(ShellCommand::ContainersSetExclusive {
                container_id,
                exclusive,
            })
        }
        other => Err(format!(
            "Unknown containers subcommand '{other}' (expected set-exclusive)"
        )),
    }
}

pub(super) fn parse_candidates_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "candidates requires a subcommand: list, delete, generate, generate-between-anchors, show, metrics, score, score-distance, score-weighted, top-k, pareto, filter, set-op, macro, template-list, template-show, template-put, template-delete, template-run"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "list" => {
            if tokens.len() > 2 {
                return Err("candidates list takes no options".to_string());
            }
            Ok(ShellCommand::CandidatesList)
        }
        "delete" => {
            if tokens.len() != 3 {
                return Err("candidates delete requires SET_NAME".to_string());
            }
            Ok(ShellCommand::CandidatesDelete {
                set_name: tokens[2].clone(),
            })
        }
        "generate" => {
            if tokens.len() < 4 {
                return Err(
                    "candidates generate requires SET_NAME SEQ_ID --length N [--step N] [--feature-kind KIND] [--feature-label-regex REGEX] [--max-distance N] [--feature-geometry MODE] [--feature-boundary MODE] [--strand-relation MODE] [--limit N]"
                        .to_string(),
                );
            }
            let set_name = tokens[2].clone();
            let seq_id = tokens[3].clone();
            let mut length_bp: Option<usize> = None;
            let mut step_bp: usize = 1;
            let mut feature_kinds: Vec<String> = vec![];
            let mut feature_label_regex: Option<String> = None;
            let mut max_distance_bp: Option<usize> = None;
            let mut feature_geometry_mode: Option<CandidateFeatureGeometryMode> = None;
            let mut feature_boundary_mode: Option<CandidateFeatureBoundaryMode> = None;
            let mut feature_strand_relation: Option<CandidateFeatureStrandRelation> = None;
            let mut limit: usize = DEFAULT_CANDIDATE_SET_LIMIT;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--length" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--length", "candidates generate")?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --length value '{raw}': {e}"))?;
                        if parsed == 0 {
                            return Err("--length must be >= 1".to_string());
                        }
                        length_bp = Some(parsed);
                    }
                    "--step" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--step", "candidates generate")?;
                        step_bp = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --step value '{raw}': {e}"))?;
                        if step_bp == 0 {
                            return Err("--step must be >= 1".to_string());
                        }
                    }
                    "--feature-kind" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-kind",
                            "candidates generate",
                        )?;
                        let trimmed = raw.trim();
                        if !trimmed.is_empty() {
                            feature_kinds.push(trimmed.to_string());
                        }
                    }
                    "--feature-label-regex" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-label-regex",
                            "candidates generate",
                        )?;
                        if raw.trim().is_empty() {
                            feature_label_regex = None;
                        } else {
                            feature_label_regex = Some(raw);
                        }
                    }
                    "--max-distance" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-distance",
                            "candidates generate",
                        )?;
                        max_distance_bp =
                            Some(raw.parse::<usize>().map_err(|e| {
                                format!("Invalid --max-distance value '{raw}': {e}")
                            })?);
                    }
                    "--feature-geometry" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-geometry",
                            "candidates generate",
                        )?;
                        feature_geometry_mode = Some(parse_candidate_feature_geometry_mode(&raw)?);
                    }
                    "--feature-boundary" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-boundary",
                            "candidates generate",
                        )?;
                        feature_boundary_mode = Some(parse_candidate_feature_boundary_mode(&raw)?);
                    }
                    "--strand-relation" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--strand-relation",
                            "candidates generate",
                        )?;
                        feature_strand_relation =
                            Some(parse_candidate_feature_strand_relation(&raw)?);
                    }
                    "--limit" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--limit", "candidates generate")?;
                        limit = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --limit value '{raw}': {e}"))?;
                        if limit == 0 {
                            return Err("--limit must be >= 1".to_string());
                        }
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for candidates generate"));
                    }
                }
            }
            let Some(length_bp) = length_bp else {
                return Err("candidates generate requires --length N".to_string());
            };
            Ok(ShellCommand::CandidatesGenerate {
                set_name,
                seq_id,
                length_bp,
                step_bp,
                feature_kinds,
                feature_label_regex,
                max_distance_bp,
                feature_geometry_mode,
                feature_boundary_mode,
                feature_strand_relation,
                limit,
            })
        }
        "generate-between-anchors" | "generate-between" => {
            if tokens.len() < 4 {
                return Err(
                    "candidates generate-between-anchors requires SET_NAME SEQ_ID --length N (--anchor-a-pos N | --anchor-a-json JSON) (--anchor-b-pos N | --anchor-b-json JSON) [--step N] [--limit N]"
                        .to_string(),
                );
            }
            let set_name = tokens[2].clone();
            let seq_id = tokens[3].clone();
            let mut length_bp: Option<usize> = None;
            let mut step_bp: usize = 1;
            let mut limit: usize = DEFAULT_CANDIDATE_SET_LIMIT;
            let mut anchor_a: Option<SequenceAnchor> = None;
            let mut anchor_b: Option<SequenceAnchor> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--length" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--length",
                            "candidates generate-between-anchors",
                        )?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --length value '{raw}': {e}"))?;
                        if parsed == 0 {
                            return Err("--length must be >= 1".to_string());
                        }
                        length_bp = Some(parsed);
                    }
                    "--step" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--step",
                            "candidates generate-between-anchors",
                        )?;
                        step_bp = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --step value '{raw}': {e}"))?;
                        if step_bp == 0 {
                            return Err("--step must be >= 1".to_string());
                        }
                    }
                    "--limit" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--limit",
                            "candidates generate-between-anchors",
                        )?;
                        limit = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --limit value '{raw}': {e}"))?;
                        if limit == 0 {
                            return Err("--limit must be >= 1".to_string());
                        }
                    }
                    "--anchor-a-pos" => {
                        if anchor_a.is_some() {
                            return Err("Anchor A was already specified".to_string());
                        }
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--anchor-a-pos",
                            "candidates generate-between-anchors",
                        )?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --anchor-a-pos value '{raw}': {e}"))?;
                        anchor_a = Some(SequenceAnchor::Position { zero_based: parsed });
                    }
                    "--anchor-b-pos" => {
                        if anchor_b.is_some() {
                            return Err("Anchor B was already specified".to_string());
                        }
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--anchor-b-pos",
                            "candidates generate-between-anchors",
                        )?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --anchor-b-pos value '{raw}': {e}"))?;
                        anchor_b = Some(SequenceAnchor::Position { zero_based: parsed });
                    }
                    "--anchor-a-json" => {
                        if anchor_a.is_some() {
                            return Err("Anchor A was already specified".to_string());
                        }
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--anchor-a-json",
                            "candidates generate-between-anchors",
                        )?;
                        anchor_a = Some(parse_sequence_anchor_json(&raw, "--anchor-a-json")?);
                    }
                    "--anchor-b-json" => {
                        if anchor_b.is_some() {
                            return Err("Anchor B was already specified".to_string());
                        }
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--anchor-b-json",
                            "candidates generate-between-anchors",
                        )?;
                        anchor_b = Some(parse_sequence_anchor_json(&raw, "--anchor-b-json")?);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for candidates generate-between-anchors"
                        ));
                    }
                }
            }
            let Some(length_bp) = length_bp else {
                return Err("candidates generate-between-anchors requires --length N".to_string());
            };
            let anchor_a = anchor_a.ok_or_else(|| {
                "candidates generate-between-anchors requires --anchor-a-pos N or --anchor-a-json JSON"
                    .to_string()
            })?;
            let anchor_b = anchor_b.ok_or_else(|| {
                "candidates generate-between-anchors requires --anchor-b-pos N or --anchor-b-json JSON"
                    .to_string()
            })?;
            Ok(ShellCommand::CandidatesGenerateBetweenAnchors {
                set_name,
                seq_id,
                anchor_a,
                anchor_b,
                length_bp,
                step_bp,
                limit,
            })
        }
        "show" => {
            if tokens.len() < 3 {
                return Err(
                    "candidates show requires SET_NAME [--limit N] [--offset N]".to_string()
                );
            }
            let set_name = tokens[2].clone();
            let mut limit = DEFAULT_CANDIDATE_PAGE_SIZE;
            let mut offset = 0usize;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--limit" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--limit", "candidates show")?;
                        limit = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --limit value '{raw}': {e}"))?;
                        if limit == 0 {
                            return Err("--limit must be >= 1".to_string());
                        }
                    }
                    "--offset" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--offset", "candidates show")?;
                        offset = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --offset value '{raw}': {e}"))?;
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for candidates show"));
                    }
                }
            }
            Ok(ShellCommand::CandidatesShow {
                set_name,
                limit,
                offset,
            })
        }
        "metrics" => {
            if tokens.len() != 3 {
                return Err("candidates metrics requires SET_NAME".to_string());
            }
            Ok(ShellCommand::CandidatesMetrics {
                set_name: tokens[2].clone(),
            })
        }
        "score" => {
            if tokens.len() < 5 {
                return Err("candidates score requires SET_NAME METRIC_NAME EXPRESSION".to_string());
            }
            let set_name = tokens[2].clone();
            let metric = tokens[3].clone();
            let expression = tokens[4..].join(" ");
            if expression.trim().is_empty() {
                return Err("candidates score requires non-empty EXPRESSION".to_string());
            }
            Ok(ShellCommand::CandidatesScoreExpression {
                set_name,
                metric,
                expression,
            })
        }
        "score-distance" => {
            if tokens.len() < 4 {
                return Err(
                    "candidates score-distance requires SET_NAME METRIC_NAME [--feature-kind KIND] [--feature-label-regex REGEX] [--feature-geometry MODE] [--feature-boundary MODE] [--strand-relation MODE]"
                        .to_string(),
                );
            }
            let set_name = tokens[2].clone();
            let metric = tokens[3].clone();
            let mut feature_kinds: Vec<String> = vec![];
            let mut feature_label_regex: Option<String> = None;
            let mut feature_geometry_mode: Option<CandidateFeatureGeometryMode> = None;
            let mut feature_boundary_mode: Option<CandidateFeatureBoundaryMode> = None;
            let mut feature_strand_relation: Option<CandidateFeatureStrandRelation> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--feature-kind" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-kind",
                            "candidates score-distance",
                        )?;
                        let trimmed = raw.trim();
                        if !trimmed.is_empty() {
                            feature_kinds.push(trimmed.to_string());
                        }
                    }
                    "--feature-label-regex" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-label-regex",
                            "candidates score-distance",
                        )?;
                        feature_label_regex = Some(raw);
                    }
                    "--feature-geometry" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-geometry",
                            "candidates score-distance",
                        )?;
                        feature_geometry_mode = Some(parse_candidate_feature_geometry_mode(&raw)?);
                    }
                    "--feature-boundary" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-boundary",
                            "candidates score-distance",
                        )?;
                        feature_boundary_mode = Some(parse_candidate_feature_boundary_mode(&raw)?);
                    }
                    "--strand-relation" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--strand-relation",
                            "candidates score-distance",
                        )?;
                        feature_strand_relation =
                            Some(parse_candidate_feature_strand_relation(&raw)?);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for candidates score-distance"
                        ));
                    }
                }
            }
            Ok(ShellCommand::CandidatesScoreDistance {
                set_name,
                metric,
                feature_kinds,
                feature_label_regex,
                feature_geometry_mode,
                feature_boundary_mode,
                feature_strand_relation,
            })
        }
        "score-weighted" => {
            if tokens.len() < 5 {
                return Err(
                    "candidates score-weighted requires SET_NAME METRIC_NAME --term METRIC:WEIGHT[:max|min] [--term ...] [--normalize|--no-normalize]"
                        .to_string(),
                );
            }
            let set_name = tokens[2].clone();
            let metric = tokens[3].clone();
            let mut objectives: Vec<CandidateWeightedObjectiveTerm> = vec![];
            let mut normalize_metrics = true;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--term" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--term",
                            "candidates score-weighted",
                        )?;
                        objectives.push(parse_weighted_objective_term(&raw)?);
                    }
                    "--normalize" => {
                        normalize_metrics = true;
                        idx += 1;
                    }
                    "--no-normalize" => {
                        normalize_metrics = false;
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for candidates score-weighted"
                        ));
                    }
                }
            }
            if objectives.is_empty() {
                return Err(
                    "candidates score-weighted requires at least one --term METRIC:WEIGHT[:max|min]"
                        .to_string(),
                );
            }
            Ok(ShellCommand::CandidatesScoreWeightedObjective {
                set_name,
                metric,
                objectives,
                normalize_metrics,
            })
        }
        "top-k" => {
            if tokens.len() < 4 {
                return Err(
                    "candidates top-k requires INPUT_SET OUTPUT_SET --metric METRIC_NAME --k N [--direction max|min] [--tie-break POLICY]"
                        .to_string(),
                );
            }
            let input_set = tokens[2].clone();
            let output_set = tokens[3].clone();
            let mut metric: Option<String> = None;
            let mut k: Option<usize> = None;
            let mut direction = CandidateObjectiveDirection::Maximize;
            let mut tie_break = CandidateTieBreakPolicy::SeqStartEnd;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--metric" => {
                        metric = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--metric",
                            "candidates top-k",
                        )?);
                    }
                    "--k" => {
                        let raw = parse_option_path(tokens, &mut idx, "--k", "candidates top-k")?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --k value '{raw}': {e}"))?;
                        if parsed == 0 {
                            return Err("--k must be >= 1".to_string());
                        }
                        k = Some(parsed);
                    }
                    "--direction" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--direction", "candidates top-k")?;
                        direction = parse_candidate_objective_direction(&raw)?;
                    }
                    "--tie-break" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--tie-break", "candidates top-k")?;
                        tie_break = parse_candidate_tie_break_policy(&raw)?;
                    }
                    other => return Err(format!("Unknown option '{other}' for candidates top-k")),
                }
            }
            let metric = metric.ok_or_else(|| "candidates top-k requires --metric".to_string())?;
            let k = k.ok_or_else(|| "candidates top-k requires --k N".to_string())?;
            Ok(ShellCommand::CandidatesTopK {
                input_set,
                output_set,
                metric,
                k,
                direction,
                tie_break,
            })
        }
        "pareto" => {
            if tokens.len() < 4 {
                return Err(
                    "candidates pareto requires INPUT_SET OUTPUT_SET --objective METRIC[:max|min] [--objective ...] [--max-candidates N] [--tie-break POLICY]"
                        .to_string(),
                );
            }
            let input_set = tokens[2].clone();
            let output_set = tokens[3].clone();
            let mut objectives: Vec<CandidateObjectiveSpec> = vec![];
            let mut max_candidates: Option<usize> = None;
            let mut tie_break = CandidateTieBreakPolicy::SeqStartEnd;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--objective" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--objective",
                            "candidates pareto",
                        )?;
                        objectives.push(parse_pareto_objective(&raw)?);
                    }
                    "--max-candidates" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-candidates",
                            "candidates pareto",
                        )?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --max-candidates value '{raw}': {e}"))?;
                        if parsed == 0 {
                            return Err("--max-candidates must be >= 1".to_string());
                        }
                        max_candidates = Some(parsed);
                    }
                    "--tie-break" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--tie-break",
                            "candidates pareto",
                        )?;
                        tie_break = parse_candidate_tie_break_policy(&raw)?;
                    }
                    other => return Err(format!("Unknown option '{other}' for candidates pareto")),
                }
            }
            if objectives.is_empty() {
                return Err(
                    "candidates pareto requires at least one --objective METRIC[:max|min]"
                        .to_string(),
                );
            }
            Ok(ShellCommand::CandidatesParetoFrontier {
                input_set,
                output_set,
                objectives,
                max_candidates,
                tie_break,
            })
        }
        "filter" => {
            if tokens.len() < 5 {
                return Err("candidates filter requires INPUT_SET OUTPUT_SET --metric METRIC_NAME [--min N] [--max N] [--min-quantile Q] [--max-quantile Q]".to_string());
            }
            let input_set = tokens[2].clone();
            let output_set = tokens[3].clone();
            let mut metric: Option<String> = None;
            let mut min: Option<f64> = None;
            let mut max: Option<f64> = None;
            let mut min_quantile: Option<f64> = None;
            let mut max_quantile: Option<f64> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--metric" => {
                        metric = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--metric",
                            "candidates filter",
                        )?);
                    }
                    "--min" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--min", "candidates filter")?;
                        min = Some(
                            raw.parse::<f64>()
                                .map_err(|e| format!("Invalid --min value '{raw}': {e}"))?,
                        );
                    }
                    "--max" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--max", "candidates filter")?;
                        max = Some(
                            raw.parse::<f64>()
                                .map_err(|e| format!("Invalid --max value '{raw}': {e}"))?,
                        );
                    }
                    "--min-quantile" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-quantile",
                            "candidates filter",
                        )?;
                        min_quantile =
                            Some(raw.parse::<f64>().map_err(|e| {
                                format!("Invalid --min-quantile value '{raw}': {e}")
                            })?);
                    }
                    "--max-quantile" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-quantile",
                            "candidates filter",
                        )?;
                        max_quantile =
                            Some(raw.parse::<f64>().map_err(|e| {
                                format!("Invalid --max-quantile value '{raw}': {e}")
                            })?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for candidates filter"));
                    }
                }
            }
            let metric = metric.ok_or_else(|| "candidates filter requires --metric".to_string())?;
            if min.zip(max).map(|(lo, hi)| lo > hi).unwrap_or(false) {
                return Err("--min must be <= --max".to_string());
            }
            if min_quantile
                .zip(max_quantile)
                .map(|(lo, hi)| lo > hi)
                .unwrap_or(false)
            {
                return Err("--min-quantile must be <= --max-quantile".to_string());
            }
            for (name, value) in [
                ("min-quantile", min_quantile),
                ("max-quantile", max_quantile),
            ] {
                if let Some(q) = value {
                    if !(0.0..=1.0).contains(&q) {
                        return Err(format!("--{name} must be between 0 and 1"));
                    }
                }
            }
            if min.is_none() && max.is_none() && min_quantile.is_none() && max_quantile.is_none() {
                return Err(
                    "candidates filter requires at least one of --min/--max/--min-quantile/--max-quantile"
                        .to_string(),
                );
            }
            Ok(ShellCommand::CandidatesFilter {
                input_set,
                output_set,
                metric,
                min,
                max,
                min_quantile,
                max_quantile,
            })
        }
        "set-op" => {
            if tokens.len() != 6 {
                return Err(
                    "candidates set-op requires union|intersect|subtract LEFT_SET RIGHT_SET OUTPUT_SET"
                        .to_string(),
                );
            }
            let op = CandidateSetOperator::parse(&tokens[2]).ok_or_else(|| {
                format!(
                    "Unsupported candidates set-op '{}'; expected union|intersect|subtract",
                    tokens[2]
                )
            })?;
            Ok(ShellCommand::CandidatesSetOp {
                op,
                left_set: tokens[3].clone(),
                right_set: tokens[4].clone(),
                output_set: tokens[5].clone(),
            })
        }
        "macro" => {
            if tokens.len() < 3 {
                return Err(
                    "candidates macro requires SCRIPT_OR_@FILE (or --file PATH), optionally with --transactional".to_string(),
                );
            }
            let mut idx = 2usize;
            let mut transactional = false;
            let mut script_file: Option<String> = None;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--transactional" | "--atomic" => {
                        transactional = true;
                        idx += 1;
                    }
                    "--file" => {
                        if script_file.is_some() {
                            return Err(
                                "candidates macro --file may only be specified once".to_string()
                            );
                        }
                        idx += 1;
                        if idx >= tokens.len() {
                            return Err("candidates macro --file requires PATH".to_string());
                        }
                        script_file = Some(tokens[idx].trim().to_string());
                        idx += 1;
                    }
                    _ => break,
                }
            }
            let script = if let Some(path) = script_file {
                if idx != tokens.len() {
                    return Err(
                        "candidates macro does not accept inline script after --file PATH"
                            .to_string(),
                    );
                }
                format!("@{path}")
            } else {
                if idx >= tokens.len() {
                    return Err("candidates macro requires SCRIPT_OR_@FILE".to_string());
                }
                tokens[idx..].join(" ")
            };
            if script.trim().is_empty() {
                return Err("candidates macro requires non-empty script".to_string());
            }
            Ok(ShellCommand::CandidatesMacro {
                script,
                transactional,
            })
        }
        "template-list" => {
            if tokens.len() != 2 {
                return Err("candidates template-list takes no options".to_string());
            }
            Ok(ShellCommand::CandidatesTemplateList)
        }
        "template-show" => {
            if tokens.len() != 3 {
                return Err("candidates template-show requires TEMPLATE_NAME".to_string());
            }
            Ok(ShellCommand::CandidatesTemplateShow {
                name: tokens[2].clone(),
            })
        }
        "template-put" | "template-upsert" => {
            if tokens.len() < 4 {
                return Err(
                    "candidates template-put requires TEMPLATE_NAME (--script SCRIPT_OR_@FILE | --file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...]"
                        .to_string(),
                );
            }
            let name = tokens[2].clone();
            let mut description: Option<String> = None;
            let mut details_url: Option<String> = None;
            let mut parameters: Vec<CandidateMacroTemplateParam> = vec![];
            let mut script: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--description" => {
                        description = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--description",
                            "candidates template-put",
                        )?);
                    }
                    "--details-url" | "--url" => {
                        if details_url.is_some() {
                            return Err(
                                "candidates template-put details URL was already specified"
                                    .to_string(),
                            );
                        }
                        details_url = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--details-url",
                            "candidates template-put",
                        )?);
                    }
                    "--param" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--param",
                            "candidates template-put",
                        )?;
                        parameters.push(parse_candidate_template_param_spec(&raw)?);
                    }
                    "--script" => {
                        if script.is_some() {
                            return Err(
                                "candidates template-put script was already specified".to_string()
                            );
                        }
                        script = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--script",
                            "candidates template-put",
                        )?);
                    }
                    "--file" => {
                        if script.is_some() {
                            return Err(
                                "candidates template-put script was already specified".to_string()
                            );
                        }
                        let path = parse_option_path(
                            tokens,
                            &mut idx,
                            "--file",
                            "candidates template-put",
                        )?;
                        script = Some(format!("@{path}"));
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for candidates template-put"
                        ));
                    }
                }
            }
            let script = script.ok_or_else(|| {
                "candidates template-put requires --script SCRIPT_OR_@FILE or --file PATH"
                    .to_string()
            })?;
            Ok(ShellCommand::CandidatesTemplateUpsert {
                name,
                description,
                details_url,
                parameters,
                script,
            })
        }
        "template-delete" => {
            if tokens.len() != 3 {
                return Err("candidates template-delete requires TEMPLATE_NAME".to_string());
            }
            Ok(ShellCommand::CandidatesTemplateDelete {
                name: tokens[2].clone(),
            })
        }
        "template-run" => {
            if tokens.len() < 3 {
                return Err(
                    "candidates template-run requires TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional]"
                        .to_string(),
                );
            }
            let name = tokens[2].clone();
            let mut bindings: HashMap<String, String> = HashMap::new();
            let mut transactional = false;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--bind" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--bind",
                            "candidates template-run",
                        )?;
                        let (key, value) = parse_template_binding(&raw)?;
                        if bindings.insert(key.clone(), value).is_some() {
                            return Err(format!(
                                "Duplicate --bind key '{}' in candidates template-run",
                                key
                            ));
                        }
                    }
                    "--transactional" | "--atomic" => {
                        transactional = true;
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for candidates template-run"
                        ));
                    }
                }
            }
            Ok(ShellCommand::CandidatesTemplateRun {
                name,
                bindings,
                transactional,
            })
        }
        other => Err(format!(
            "Unknown candidates subcommand '{other}' (expected list, delete, generate, generate-between-anchors, show, metrics, score, score-distance, score-weighted, top-k, pareto, filter, set-op, macro, template-list, template-show, template-put, template-delete, template-run)"
        )),
    }
}

pub(super) fn parse_guides_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "guides requires a subcommand: list, show, put, delete, filter, filter-show, oligos-generate, oligos-list, oligos-show, oligos-export, protocol-export"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "list" => {
            if tokens.len() > 2 {
                return Err("guides list takes no options".to_string());
            }
            Ok(ShellCommand::GuidesList)
        }
        "show" => {
            if tokens.len() < 3 {
                return Err(
                    "guides show requires GUIDE_SET_ID [--limit N] [--offset N]".to_string()
                );
            }
            let guide_set_id = tokens[2].clone();
            let mut limit = DEFAULT_GUIDE_PAGE_SIZE;
            let mut offset = 0usize;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--limit" => {
                        let raw = parse_option_path(tokens, &mut idx, "--limit", "guides show")?;
                        limit = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --limit value '{raw}': {e}"))?;
                        if limit == 0 {
                            return Err("--limit must be >= 1".to_string());
                        }
                    }
                    "--offset" => {
                        let raw = parse_option_path(tokens, &mut idx, "--offset", "guides show")?;
                        offset = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --offset value '{raw}': {e}"))?;
                    }
                    other => return Err(format!("Unknown option '{other}' for guides show")),
                }
            }
            Ok(ShellCommand::GuidesShow {
                guide_set_id,
                limit,
                offset,
            })
        }
        "put" | "upsert" => {
            if tokens.len() < 4 {
                return Err(
                    "guides put requires GUIDE_SET_ID (--json JSON|@FILE | --file PATH)"
                        .to_string(),
                );
            }
            let guide_set_id = tokens[2].clone();
            let mut guides_json: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--json" => {
                        if guides_json.is_some() {
                            return Err("guides put JSON payload was already specified".to_string());
                        }
                        guides_json =
                            Some(parse_option_path(tokens, &mut idx, "--json", "guides put")?);
                    }
                    "--file" => {
                        if guides_json.is_some() {
                            return Err("guides put JSON payload was already specified".to_string());
                        }
                        let path = parse_option_path(tokens, &mut idx, "--file", "guides put")?;
                        guides_json = Some(format!("@{path}"));
                    }
                    other => return Err(format!("Unknown option '{other}' for guides put")),
                }
            }
            let guides_json = guides_json.ok_or_else(|| {
                "guides put requires --json JSON|@FILE or --file PATH".to_string()
            })?;
            Ok(ShellCommand::GuidesPut {
                guide_set_id,
                guides_json,
            })
        }
        "delete" => {
            if tokens.len() != 3 {
                return Err("guides delete requires GUIDE_SET_ID".to_string());
            }
            Ok(ShellCommand::GuidesDelete {
                guide_set_id: tokens[2].clone(),
            })
        }
        "filter" => {
            if tokens.len() < 3 {
                return Err(
                    "guides filter requires GUIDE_SET_ID [--config JSON|@FILE] [--config-file PATH] [--output-set GUIDE_SET_ID]"
                        .to_string(),
                );
            }
            let guide_set_id = tokens[2].clone();
            let mut config_json: Option<String> = None;
            let mut output_guide_set_id: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--config" => {
                        if config_json.is_some() {
                            return Err("guides filter config was already specified".to_string());
                        }
                        config_json = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--config",
                            "guides filter",
                        )?);
                    }
                    "--config-file" => {
                        if config_json.is_some() {
                            return Err("guides filter config was already specified".to_string());
                        }
                        let path =
                            parse_option_path(tokens, &mut idx, "--config-file", "guides filter")?;
                        config_json = Some(format!("@{path}"));
                    }
                    "--output-set" => {
                        output_guide_set_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--output-set",
                            "guides filter",
                        )?);
                    }
                    other => return Err(format!("Unknown option '{other}' for guides filter")),
                }
            }
            Ok(ShellCommand::GuidesFilter {
                guide_set_id,
                config_json,
                output_guide_set_id,
            })
        }
        "filter-show" => {
            if tokens.len() != 3 {
                return Err("guides filter-show requires GUIDE_SET_ID".to_string());
            }
            Ok(ShellCommand::GuidesFilterShow {
                guide_set_id: tokens[2].clone(),
            })
        }
        "oligos-generate" => {
            if tokens.len() < 4 {
                return Err(
                    "guides oligos-generate requires GUIDE_SET_ID TEMPLATE_ID [--apply-5prime-g-extension] [--output-oligo-set ID] [--passed-only]"
                        .to_string(),
                );
            }
            let guide_set_id = tokens[2].clone();
            let template_id = tokens[3].clone();
            let mut apply_5prime_g_extension = false;
            let mut output_oligo_set_id: Option<String> = None;
            let mut passed_only = false;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--apply-5prime-g-extension" => {
                        apply_5prime_g_extension = true;
                        idx += 1;
                    }
                    "--output-oligo-set" => {
                        output_oligo_set_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--output-oligo-set",
                            "guides oligos-generate",
                        )?);
                    }
                    "--passed-only" => {
                        passed_only = true;
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for guides oligos-generate"
                        ));
                    }
                }
            }
            Ok(ShellCommand::GuidesOligosGenerate {
                guide_set_id,
                template_id,
                apply_5prime_g_extension,
                output_oligo_set_id,
                passed_only,
            })
        }
        "oligos-list" => {
            let mut guide_set_id: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--guide-set" => {
                        guide_set_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--guide-set",
                            "guides oligos-list",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for guides oligos-list"));
                    }
                }
            }
            Ok(ShellCommand::GuidesOligosList { guide_set_id })
        }
        "oligos-show" => {
            if tokens.len() != 3 {
                return Err("guides oligos-show requires OLIGO_SET_ID".to_string());
            }
            Ok(ShellCommand::GuidesOligosShow {
                oligo_set_id: tokens[2].clone(),
            })
        }
        "oligos-export" => {
            if tokens.len() < 4 {
                return Err(
                    "guides oligos-export requires GUIDE_SET_ID OUTPUT_PATH [--format csv_table|plate_csv|fasta] [--plate 96|384] [--oligo-set ID]"
                        .to_string(),
                );
            }
            let guide_set_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut format = GuideOligoExportFormat::CsvTable;
            let mut plate_format: Option<GuideOligoPlateFormat> = None;
            let mut oligo_set_id: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--format" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--format",
                            "guides oligos-export",
                        )?;
                        format = parse_guide_export_format(&raw)?;
                    }
                    "--plate" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--plate", "guides oligos-export")?;
                        plate_format = Some(parse_guide_plate_format(&raw)?);
                    }
                    "--oligo-set" => {
                        oligo_set_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--oligo-set",
                            "guides oligos-export",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for guides oligos-export"));
                    }
                }
            }
            Ok(ShellCommand::GuidesOligosExport {
                guide_set_id,
                oligo_set_id,
                format,
                path,
                plate_format,
            })
        }
        "protocol-export" => {
            if tokens.len() < 4 {
                return Err(
                    "guides protocol-export requires GUIDE_SET_ID OUTPUT_PATH [--oligo-set ID] [--no-qc]"
                        .to_string(),
                );
            }
            let guide_set_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut oligo_set_id: Option<String> = None;
            let mut include_qc_checklist = true;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--oligo-set" => {
                        oligo_set_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--oligo-set",
                            "guides protocol-export",
                        )?);
                    }
                    "--no-qc" => {
                        include_qc_checklist = false;
                        idx += 1;
                    }
                    "--with-qc" => {
                        include_qc_checklist = true;
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for guides protocol-export"
                        ));
                    }
                }
            }
            Ok(ShellCommand::GuidesProtocolExport {
                guide_set_id,
                oligo_set_id,
                path,
                include_qc_checklist,
            })
        }
        other => Err(format!(
            "Unknown guides subcommand '{other}' (expected list, show, put, delete, filter, filter-show, oligos-generate, oligos-list, oligos-show, oligos-export, protocol-export)"
        )),
    }
}

pub(super) fn sequence_feature_roi_range_0based(
    dna: &DNAsequence,
    feature_id: usize,
) -> Result<(usize, usize), String> {
    let feature = dna
        .features()
        .get(feature_id)
        .ok_or_else(|| format!("Feature id {feature_id} is out of range"))?;
    let mut ranges: Vec<(usize, usize)> = Vec::new();
    collect_location_ranges_usize(&feature.location, &mut ranges);
    if ranges.is_empty() {
        let (from, to) = feature
            .location
            .find_bounds()
            .map_err(|e| format!("Could not read bounds for feature id {feature_id}: {e}"))?;
        if from < 0 || to < 0 {
            return Err(format!("Feature id {feature_id} has negative bounds"));
        }
        ranges.push((from as usize, to as usize));
    }
    let seq_len = dna.len();
    if seq_len == 0 {
        return Err("Template sequence is empty".to_string());
    }
    let start = ranges
        .iter()
        .map(|(start, _)| *start)
        .min()
        .ok_or_else(|| format!("Feature id {feature_id} has no usable ranges"))?;
    if start >= seq_len {
        return Err(format!(
            "Feature id {feature_id} starts outside sequence length {seq_len}"
        ));
    }
    let end_exclusive = ranges
        .iter()
        .map(|(_, end)| *end)
        .max()
        .ok_or_else(|| format!("Feature id {feature_id} has no usable ranges"))?
        .min(seq_len);
    if end_exclusive <= start {
        return Err(format!(
            "Feature id {feature_id} has invalid range {}..{} for sequence length {seq_len}",
            start, end_exclusive
        ));
    }
    Ok((start, end_exclusive))
}

fn parse_feature_range(raw: &str, context: &str) -> Result<(usize, usize), String> {
    let trimmed = raw.trim();
    let (left, right) = if let Some((left, right)) = trimmed.split_once("..") {
        (left.trim(), right.trim())
    } else if let Some((left, right)) = trimmed.split_once(':') {
        (left.trim(), right.trim())
    } else {
        return Err(format!(
            "Invalid --range value '{raw}' for {context}; expected START..END"
        ));
    };
    let start = left
        .parse::<usize>()
        .map_err(|e| format!("Invalid range start '{left}' for {context}: {e}"))?;
    let end_exclusive = right
        .parse::<usize>()
        .map_err(|e| format!("Invalid range end '{right}' for {context}: {e}"))?;
    if end_exclusive <= start {
        return Err(format!(
            "Invalid --range value '{raw}' for {context}: end must be > start"
        ));
    }
    Ok((start, end_exclusive))
}

fn parse_feature_query_strand(raw: &str) -> Result<SequenceFeatureStrandFilter, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "any" => Ok(SequenceFeatureStrandFilter::Any),
        "forward" | "+" | "plus" => Ok(SequenceFeatureStrandFilter::Forward),
        "reverse" | "-" | "minus" => Ok(SequenceFeatureStrandFilter::Reverse),
        other => Err(format!(
            "Unsupported --strand value '{other}' (expected any|forward|reverse)"
        )),
    }
}

fn parse_feature_query_sort(raw: &str) -> Result<SequenceFeatureSortBy, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "feature_id" | "id" => Ok(SequenceFeatureSortBy::FeatureId),
        "start" => Ok(SequenceFeatureSortBy::Start),
        "end" => Ok(SequenceFeatureSortBy::End),
        "kind" => Ok(SequenceFeatureSortBy::Kind),
        "length" | "len" => Ok(SequenceFeatureSortBy::Length),
        other => Err(format!(
            "Unsupported --sort value '{other}' (expected feature_id|start|end|kind|length)"
        )),
    }
}

fn parse_feature_qual_key_value(raw: &str, flag: &str) -> Result<(String, String), String> {
    let Some((key, value)) = raw.split_once('=') else {
        return Err(format!("{flag} requires KEY=VALUE (received '{raw}')"));
    };
    let key = key.trim().to_string();
    let value = value.trim().to_string();
    if key.is_empty() {
        return Err(format!("{flag} requires non-empty KEY"));
    }
    if value.is_empty() {
        return Err(format!("{flag} requires non-empty VALUE"));
    }
    Ok((key, value))
}

#[derive(Default)]
struct FeatureQueryOptionState {
    range_arg: Option<(usize, usize)>,
    start_arg: Option<usize>,
    end_arg: Option<usize>,
}

fn parse_feature_bed_coordinate_mode(raw: &str) -> Result<FeatureBedCoordinateMode, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "auto" => Ok(FeatureBedCoordinateMode::Auto),
        "local" => Ok(FeatureBedCoordinateMode::Local),
        "genomic" => Ok(FeatureBedCoordinateMode::Genomic),
        other => Err(format!(
            "Unsupported --coordinate-mode value '{other}' (expected auto|local|genomic)"
        )),
    }
}

fn parse_inline_sequence_topology(raw: &str) -> Result<InlineSequenceTopology, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "linear" | "lin" => Ok(InlineSequenceTopology::Linear),
        "circular" | "circ" => Ok(InlineSequenceTopology::Circular),
        other => Err(format!(
            "Unsupported --topology value '{other}' (expected linear|circular)"
        )),
    }
}

fn parse_tf_threshold_override_value(
    raw: &str,
    flag: &str,
    context: &str,
) -> Result<(String, f64), String> {
    let Some((tf, value)) = raw.split_once('=') else {
        return Err(format!(
            "{flag} for {context} requires TF=VALUE (received '{raw}')"
        ));
    };
    let tf = tf.trim().to_string();
    if tf.is_empty() {
        return Err(format!("{flag} for {context} requires non-empty TF"));
    }
    let value = value
        .trim()
        .parse::<f64>()
        .map_err(|e| format!("Invalid {flag} value '{raw}' for {context}: {e}"))?;
    Ok((tf, value))
}

fn try_parse_feature_query_option(
    tokens: &[String],
    idx: &mut usize,
    context: &str,
    query: &mut SequenceFeatureQuery,
    state: &mut FeatureQueryOptionState,
) -> Result<bool, String> {
    match tokens[*idx].as_str() {
        "--kind" => {
            let raw = parse_option_path(tokens, idx, "--kind", context)?;
            let trimmed = raw.trim();
            if !trimmed.is_empty() {
                query.kind_in.push(trimmed.to_string());
            }
            Ok(true)
        }
        "--kind-not" => {
            let raw = parse_option_path(tokens, idx, "--kind-not", context)?;
            let trimmed = raw.trim();
            if !trimmed.is_empty() {
                query.kind_not_in.push(trimmed.to_string());
            }
            Ok(true)
        }
        "--range" => {
            let raw = parse_option_path(tokens, idx, "--range", context)?;
            if state.range_arg.is_some() {
                return Err("--range was specified multiple times".to_string());
            }
            state.range_arg = Some(parse_feature_range(&raw, context)?);
            Ok(true)
        }
        "--start" => {
            let raw = parse_option_path(tokens, idx, "--start", context)?;
            let parsed = raw
                .parse::<usize>()
                .map_err(|e| format!("Invalid --start value '{raw}': {e}"))?;
            state.start_arg = Some(parsed);
            Ok(true)
        }
        "--end" => {
            let raw = parse_option_path(tokens, idx, "--end", context)?;
            let parsed = raw
                .parse::<usize>()
                .map_err(|e| format!("Invalid --end value '{raw}': {e}"))?;
            state.end_arg = Some(parsed);
            Ok(true)
        }
        "--overlap" => {
            query.range_relation = SequenceFeatureRangeRelation::Overlap;
            *idx += 1;
            Ok(true)
        }
        "--within" => {
            query.range_relation = SequenceFeatureRangeRelation::Within;
            *idx += 1;
            Ok(true)
        }
        "--contains" => {
            query.range_relation = SequenceFeatureRangeRelation::Contains;
            *idx += 1;
            Ok(true)
        }
        "--strand" => {
            let raw = parse_option_path(tokens, idx, "--strand", context)?;
            query.strand = parse_feature_query_strand(&raw)?;
            Ok(true)
        }
        "--label" => {
            query.label_contains = Some(parse_option_path(tokens, idx, "--label", context)?);
            Ok(true)
        }
        "--label-regex" => {
            query.label_regex = Some(parse_option_path(tokens, idx, "--label-regex", context)?);
            Ok(true)
        }
        "--qual" => {
            let key = parse_option_path(tokens, idx, "--qual", context)?;
            query
                .qualifier_filters
                .push(SequenceFeatureQualifierFilter {
                    key,
                    ..SequenceFeatureQualifierFilter::default()
                });
            Ok(true)
        }
        "--qual-contains" => {
            let raw = parse_option_path(tokens, idx, "--qual-contains", context)?;
            let (key, value_contains) = parse_feature_qual_key_value(&raw, "--qual-contains")?;
            query
                .qualifier_filters
                .push(SequenceFeatureQualifierFilter {
                    key,
                    value_contains: Some(value_contains),
                    ..SequenceFeatureQualifierFilter::default()
                });
            Ok(true)
        }
        "--qual-regex" => {
            let raw = parse_option_path(tokens, idx, "--qual-regex", context)?;
            let (key, value_regex) = parse_feature_qual_key_value(&raw, "--qual-regex")?;
            query
                .qualifier_filters
                .push(SequenceFeatureQualifierFilter {
                    key,
                    value_regex: Some(value_regex),
                    ..SequenceFeatureQualifierFilter::default()
                });
            Ok(true)
        }
        "--min-len" => {
            let raw = parse_option_path(tokens, idx, "--min-len", context)?;
            let parsed = raw
                .parse::<usize>()
                .map_err(|e| format!("Invalid --min-len value '{raw}': {e}"))?;
            query.min_len_bp = Some(parsed);
            Ok(true)
        }
        "--max-len" => {
            let raw = parse_option_path(tokens, idx, "--max-len", context)?;
            let parsed = raw
                .parse::<usize>()
                .map_err(|e| format!("Invalid --max-len value '{raw}': {e}"))?;
            query.max_len_bp = Some(parsed);
            Ok(true)
        }
        "--limit" => {
            let raw = parse_option_path(tokens, idx, "--limit", context)?;
            let parsed = raw
                .parse::<usize>()
                .map_err(|e| format!("Invalid --limit value '{raw}': {e}"))?;
            if parsed == 0 {
                return Err("--limit must be >= 1".to_string());
            }
            query.limit = Some(parsed);
            Ok(true)
        }
        "--offset" => {
            let raw = parse_option_path(tokens, idx, "--offset", context)?;
            query.offset = raw
                .parse::<usize>()
                .map_err(|e| format!("Invalid --offset value '{raw}': {e}"))?;
            Ok(true)
        }
        "--sort" => {
            let raw = parse_option_path(tokens, idx, "--sort", context)?;
            query.sort_by = parse_feature_query_sort(&raw)?;
            Ok(true)
        }
        "--desc" => {
            query.descending = true;
            *idx += 1;
            Ok(true)
        }
        "--asc" => {
            query.descending = false;
            *idx += 1;
            Ok(true)
        }
        "--include-source" => {
            query.include_source = true;
            *idx += 1;
            Ok(true)
        }
        "--include-qualifiers" => {
            query.include_qualifiers = true;
            *idx += 1;
            Ok(true)
        }
        _ => Ok(false),
    }
}

fn finalize_feature_query_options(
    state: FeatureQueryOptionState,
    context: &str,
    query: &mut SequenceFeatureQuery,
) -> Result<(), String> {
    if state.range_arg.is_some() && (state.start_arg.is_some() || state.end_arg.is_some()) {
        return Err(format!(
            "{context} accepts either --range START..END or --start/--end, not both"
        ));
    }
    if let Some((start, end_exclusive)) = state.range_arg {
        query.start_0based = Some(start);
        query.end_0based_exclusive = Some(end_exclusive);
    } else if state.start_arg.is_some() || state.end_arg.is_some() {
        let start = state
            .start_arg
            .ok_or_else(|| format!("{context} --start requires matching --end"))?;
        let end_exclusive = state
            .end_arg
            .ok_or_else(|| format!("{context} --end requires matching --start"))?;
        if end_exclusive <= start {
            return Err("--end must be > --start".to_string());
        }
        query.start_0based = Some(start);
        query.end_0based_exclusive = Some(end_exclusive);
    }
    Ok(())
}

fn build_sequence_scan_target_from_feature_state(
    seq_id: Option<String>,
    sequence_text: Option<String>,
    topology: InlineSequenceTopology,
    id_hint: Option<String>,
    state: FeatureQueryOptionState,
    context: &str,
) -> Result<SequenceScanTarget, String> {
    let (span_start_0based, span_end_0based_exclusive) =
        if state.range_arg.is_some() || state.start_arg.is_some() || state.end_arg.is_some() {
            let mut query = SequenceFeatureQuery::default();
            finalize_feature_query_options(state, context, &mut query)?;
            (query.start_0based, query.end_0based_exclusive)
        } else {
            (None, None)
        };
    match (seq_id, sequence_text) {
        (Some(seq_id), None) => Ok(SequenceScanTarget::SeqId {
            seq_id,
            span_start_0based,
            span_end_0based_exclusive,
        }),
        (None, Some(sequence_text)) => Ok(SequenceScanTarget::InlineSequence {
            sequence_text,
            topology,
            id_hint,
            span_start_0based,
            span_end_0based_exclusive,
        }),
        (Some(_), Some(_)) => Err(format!(
            "{context} accepts either SEQ_ID or --sequence-text DNA, not both"
        )),
        (None, None) => Err(format!(
            "{context} requires either SEQ_ID or --sequence-text DNA"
        )),
    }
}

pub(super) fn parse_features_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "features requires a subcommand: query, export-bed, tfbs-summary, tfbs-score-tracks-svg, tfbs-track-similarity, tfbs-score-track-correlation-svg, tfbs-scan, restriction-scan"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "query" => {
            if tokens.len() < 3 {
                return Err(
                    "features query requires SEQ_ID [--kind KIND] [--kind-not KIND] [--range START..END|--start N --end N] [--overlap|--within|--contains] [--strand any|forward|reverse] [--label TEXT] [--label-regex REGEX] [--qual KEY] [--qual-contains KEY=VALUE] [--qual-regex KEY=REGEX] [--min-len N] [--max-len N] [--limit N] [--offset N] [--sort feature_id|start|end|kind|length] [--desc] [--include-source] [--include-qualifiers]"
                        .to_string(),
                );
            }
            let mut query = SequenceFeatureQuery {
                seq_id: tokens[2].clone(),
                ..SequenceFeatureQuery::default()
            };
            let mut state = FeatureQueryOptionState::default();
            let mut idx = 3usize;
            while idx < tokens.len() {
                if !try_parse_feature_query_option(
                    tokens,
                    &mut idx,
                    "features query",
                    &mut query,
                    &mut state,
                )? {
                    return Err(format!(
                        "Unknown option '{}' for features query",
                        tokens[idx]
                    ));
                }
            }
            finalize_feature_query_options(state, "features query", &mut query)?;
            Ok(ShellCommand::FeaturesQuery { query })
        }
        "export-bed" => {
            if tokens.len() < 4 {
                return Err(
                    "features export-bed requires SEQ_ID OUTPUT.bed [--coordinate-mode auto|local|genomic] [--include-restriction-sites] [--restriction-enzyme NAME] [--kind KIND] [--kind-not KIND] [--range START..END|--start N --end N] [--overlap|--within|--contains] [--strand any|forward|reverse] [--label TEXT] [--label-regex REGEX] [--qual KEY] [--qual-contains KEY=VALUE] [--qual-regex KEY=REGEX] [--min-len N] [--max-len N] [--limit N] [--offset N] [--sort feature_id|start|end|kind|length] [--desc] [--include-source] [--include-qualifiers]"
                        .to_string(),
                );
            }
            let mut query = SequenceFeatureQuery {
                seq_id: tokens[2].clone(),
                ..SequenceFeatureQuery::default()
            };
            let output = tokens[3].clone();
            let mut coordinate_mode = FeatureBedCoordinateMode::Auto;
            let mut include_restriction_sites = false;
            let mut restriction_enzymes: Vec<String> = vec![];
            let mut state = FeatureQueryOptionState::default();
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--coordinate-mode" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--coordinate-mode",
                            "features export-bed",
                        )?;
                        coordinate_mode = parse_feature_bed_coordinate_mode(&raw)?;
                    }
                    "--include-restriction-sites" => {
                        include_restriction_sites = true;
                        idx += 1;
                    }
                    "--restriction-enzyme" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--restriction-enzyme",
                            "features export-bed",
                        )?;
                        let trimmed = raw.trim();
                        if !trimmed.is_empty() {
                            restriction_enzymes.push(trimmed.to_string());
                        }
                    }
                    other => {
                        if !try_parse_feature_query_option(
                            tokens,
                            &mut idx,
                            "features export-bed",
                            &mut query,
                            &mut state,
                        )? {
                            return Err(format!(
                                "Unknown option '{other}' for features export-bed"
                            ));
                        }
                    }
                }
            }
            finalize_feature_query_options(state, "features export-bed", &mut query)?;
            Ok(ShellCommand::FeaturesExportBed {
                query,
                output,
                coordinate_mode,
                include_restriction_sites,
                restriction_enzymes,
            })
        }
        "tfbs-summary" => {
            if tokens.len() < 3 {
                return Err(
                    "features tfbs-summary requires SEQ_ID --focus START..END [--context START..END] [--min-focus-count N] [--min-context-count N] [--limit N]"
                        .to_string(),
                );
            }
            let mut request = TfbsRegionSummaryRequest {
                seq_id: tokens[2].clone(),
                ..TfbsRegionSummaryRequest::default()
            };
            let mut focus_range: Option<(usize, usize)> = None;
            let mut context_range: Option<(usize, usize)> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--focus" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--focus",
                            "features tfbs-summary",
                        )?;
                        if focus_range.is_some() {
                            return Err("--focus was specified multiple times".to_string());
                        }
                        focus_range = Some(parse_feature_range(&raw, "features tfbs-summary")?);
                    }
                    "--context" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--context",
                            "features tfbs-summary",
                        )?;
                        if context_range.is_some() {
                            return Err("--context was specified multiple times".to_string());
                        }
                        context_range = Some(parse_feature_range(&raw, "features tfbs-summary")?);
                    }
                    "--min-focus-count" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-focus-count",
                            "features tfbs-summary",
                        )?;
                        request.min_focus_occurrences = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --min-focus-count value '{raw}': {e}"))?;
                    }
                    "--min-context-count" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-context-count",
                            "features tfbs-summary",
                        )?;
                        request.min_context_occurrences = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --min-context-count value '{raw}': {e}")
                        })?;
                    }
                    "--limit" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--limit",
                            "features tfbs-summary",
                        )?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --limit value '{raw}': {e}"))?;
                        if parsed == 0 {
                            return Err("--limit must be >= 1".to_string());
                        }
                        request.limit = Some(parsed);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for features tfbs-summary"
                        ));
                    }
                }
            }
            let (focus_start, focus_end) = focus_range
                .ok_or_else(|| "features tfbs-summary requires --focus START..END".to_string())?;
            request.focus_start_0based = focus_start;
            request.focus_end_0based_exclusive = focus_end;
            if let Some((context_start, context_end)) = context_range {
                request.context_start_0based = Some(context_start);
                request.context_end_0based_exclusive = Some(context_end);
            }
            Ok(ShellCommand::FeaturesTfbsSummary { request })
        }
        "tfbs-score-tracks-svg" => {
            if tokens.len() < 4 {
                return Err(
                    "features tfbs-score-tracks-svg requires either SEQ_ID OUTPUT.svg or --sequence-text DNA --output OUTPUT.svg [--topology linear|circular] [--id-hint TEXT] --motif TOKEN [--motif TOKEN ...] [--range START..END|--start N --end N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative]"
                        .to_string(),
                );
            }
            let mut seq_id: Option<String> = None;
            let mut sequence_text: Option<String> = None;
            let mut topology = InlineSequenceTopology::Linear;
            let mut id_hint: Option<String> = None;
            let mut output: Option<String> = None;
            let mut motifs: Vec<String> = vec![];
            let mut score_kind = TfbsScoreTrackValueKind::LlrBits;
            let mut clip_negative = true;
            let mut state = FeatureQueryOptionState::default();
            let mut idx = 2usize;
            if tokens[2].starts_with("--") {
                // options-only inline form
            } else {
                let raw = tokens[2].trim().to_string();
                if raw.is_empty() {
                    return Err(
                        "features tfbs-score-tracks-svg requires non-empty SEQ_ID".to_string()
                    );
                }
                seq_id = Some(raw);
                output = Some(tokens[3].clone());
                idx = 4usize;
            }
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--sequence-text" => {
                        sequence_text = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--sequence-text",
                            "features tfbs-score-tracks-svg",
                        )?);
                    }
                    "--topology" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--topology",
                            "features tfbs-score-tracks-svg",
                        )?;
                        topology = parse_inline_sequence_topology(&raw)?;
                    }
                    "--id-hint" => {
                        id_hint = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--id-hint",
                            "features tfbs-score-tracks-svg",
                        )?);
                    }
                    "--output" => {
                        let value = parse_option_path(
                            tokens,
                            &mut idx,
                            "--output",
                            "features tfbs-score-tracks-svg",
                        )?;
                        if output.replace(value).is_some() {
                            return Err("--output was specified multiple times".to_string());
                        }
                    }
                    "--motif" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--motif",
                            "features tfbs-score-tracks-svg",
                        )?;
                        let normalized = raw.trim();
                        if normalized.is_empty() {
                            return Err("--motif must not be empty".to_string());
                        }
                        motifs.push(normalized.to_string());
                    }
                    "--motifs" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--motifs",
                            "features tfbs-score-tracks-svg",
                        )?;
                        for token in raw
                            .split(',')
                            .map(str::trim)
                            .filter(|value| !value.is_empty())
                        {
                            motifs.push(token.to_string());
                        }
                    }
                    "--range" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--range",
                            "features tfbs-score-tracks-svg",
                        )?;
                        if state.range_arg.is_some() {
                            return Err("--range was specified multiple times".to_string());
                        }
                        state.range_arg =
                            Some(parse_feature_range(&raw, "features tfbs-score-tracks-svg")?);
                    }
                    "--start" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--start",
                            "features tfbs-score-tracks-svg",
                        )?;
                        state.start_arg = Some(
                            raw.parse::<usize>()
                                .map_err(|e| format!("Invalid --start value '{raw}': {e}"))?,
                        );
                    }
                    "--end" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--end",
                            "features tfbs-score-tracks-svg",
                        )?;
                        state.end_arg = Some(
                            raw.parse::<usize>()
                                .map_err(|e| format!("Invalid --end value '{raw}': {e}"))?,
                        );
                    }
                    "--score-kind" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--score-kind",
                            "features tfbs-score-tracks-svg",
                        )?;
                        score_kind = match raw.trim() {
                            "llr_bits" => TfbsScoreTrackValueKind::LlrBits,
                            "llr_quantile" => TfbsScoreTrackValueKind::LlrQuantile,
                            "llr_background_quantile" => {
                                TfbsScoreTrackValueKind::LlrBackgroundQuantile
                            }
                            "llr_background_tail_log10" => {
                                TfbsScoreTrackValueKind::LlrBackgroundTailLog10
                            }
                            "true_log_odds_bits" => TfbsScoreTrackValueKind::TrueLogOddsBits,
                            "true_log_odds_quantile" => {
                                TfbsScoreTrackValueKind::TrueLogOddsQuantile
                            }
                            "true_log_odds_background_quantile" => {
                                TfbsScoreTrackValueKind::TrueLogOddsBackgroundQuantile
                            }
                            "true_log_odds_background_tail_log10" => {
                                TfbsScoreTrackValueKind::TrueLogOddsBackgroundTailLog10
                            }
                            other => {
                                return Err(format!(
                                    "Unsupported --score-kind value '{other}' (expected llr_bits, llr_quantile, llr_background_quantile, llr_background_tail_log10, true_log_odds_bits, true_log_odds_quantile, true_log_odds_background_quantile, or true_log_odds_background_tail_log10)"
                                ));
                            }
                        };
                    }
                    "--allow-negative" => {
                        clip_negative = false;
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for features tfbs-score-tracks-svg"
                        ));
                    }
                }
            }
            if motifs.is_empty() {
                return Err(
                    "features tfbs-score-tracks-svg requires at least one --motif TOKEN"
                        .to_string(),
                );
            }
            let output = output.ok_or_else(|| {
                "features tfbs-score-tracks-svg requires OUTPUT.svg or --output OUTPUT.svg"
                    .to_string()
            })?;
            let target = build_sequence_scan_target_from_feature_state(
                seq_id,
                sequence_text,
                topology,
                id_hint,
                state,
                "features tfbs-score-tracks-svg",
            )?;
            if !matches!(
                &target,
                SequenceScanTarget::SeqId {
                    span_end_0based_exclusive: Some(_),
                    ..
                } | SequenceScanTarget::InlineSequence {
                    span_end_0based_exclusive: Some(_),
                    ..
                }
            ) {
                return Err(
                    "features tfbs-score-tracks-svg requires --range START..END or --end N"
                        .to_string(),
                );
            }
            Ok(ShellCommand::FeaturesTfbsScoreTracksSvg {
                target,
                motifs,
                score_kind,
                clip_negative,
                output,
            })
        }
        "tfbs-track-similarity" => {
            if tokens.len() < 3 {
                return Err(
                    "features tfbs-track-similarity requires either SEQ_ID or --sequence-text DNA [--topology linear|circular] [--id-hint TEXT] --anchor-motif TOKEN [--candidate-motif TOKEN ...|--candidate-motifs CSV|--candidate-motif ALL] [--range START..END|--start N --end N] [--ranking-metric raw_pearson|smoothed_pearson|raw_spearman|smoothed_spearman] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--species TEXT] [--include-remote-metadata] [--limit N] [--path FILE.json]"
                        .to_string(),
                );
            }
            let mut seq_id: Option<String> = None;
            let mut sequence_text: Option<String> = None;
            let mut topology = InlineSequenceTopology::Linear;
            let mut id_hint: Option<String> = None;
            let mut anchor_motif: Option<String> = None;
            let mut candidate_motifs: Vec<String> = vec![];
            let mut ranking_metric = TfbsTrackSimilarityRankingMetric::SmoothedSpearman;
            let mut score_kind = TfbsScoreTrackValueKind::LlrBits;
            let mut clip_negative = true;
            let mut species_filters: Vec<String> = vec![];
            let mut include_remote_metadata = false;
            let mut limit: Option<usize> = None;
            let mut path: Option<String> = None;
            let mut state = FeatureQueryOptionState::default();
            let mut idx = 2usize;
            if tokens[2].starts_with("--") {
                // options-only inline form
            } else {
                let raw = tokens[2].trim().to_string();
                if !raw.is_empty() {
                    seq_id = Some(raw);
                }
                idx = 3usize;
            }
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--sequence-text" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--sequence-text",
                            "features tfbs-track-similarity",
                        )?;
                        if raw.trim().is_empty() {
                            return Err("--sequence-text must not be empty".to_string());
                        }
                        sequence_text = Some(raw);
                    }
                    "--topology" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--topology",
                            "features tfbs-track-similarity",
                        )?;
                        topology = match raw.trim() {
                            "linear" => InlineSequenceTopology::Linear,
                            "circular" => InlineSequenceTopology::Circular,
                            other => {
                                return Err(format!(
                                    "Unsupported --topology value '{other}' (expected linear or circular)"
                                ));
                            }
                        };
                    }
                    "--id-hint" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--id-hint",
                            "features tfbs-track-similarity",
                        )?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            id_hint = None;
                        } else {
                            id_hint = Some(trimmed.to_string());
                        }
                    }
                    "--anchor-motif" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--anchor-motif",
                            "features tfbs-track-similarity",
                        )?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            return Err("--anchor-motif must not be empty".to_string());
                        }
                        anchor_motif = Some(trimmed.to_string());
                    }
                    "--candidate-motif" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--candidate-motif",
                            "features tfbs-track-similarity",
                        )?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            return Err("--candidate-motif must not be empty".to_string());
                        }
                        candidate_motifs.push(trimmed.to_string());
                    }
                    "--candidate-motifs" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--candidate-motifs",
                            "features tfbs-track-similarity",
                        )?;
                        for token in raw
                            .split(',')
                            .map(str::trim)
                            .filter(|value| !value.is_empty())
                        {
                            candidate_motifs.push(token.to_string());
                        }
                    }
                    "--range" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--range",
                            "features tfbs-track-similarity",
                        )?;
                        if state.range_arg.is_some() {
                            return Err("--range was specified multiple times".to_string());
                        }
                        state.range_arg =
                            Some(parse_feature_range(&raw, "features tfbs-track-similarity")?);
                    }
                    "--start" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--start",
                            "features tfbs-track-similarity",
                        )?;
                        state.start_arg = Some(
                            raw.parse::<usize>()
                                .map_err(|e| format!("Invalid --start value '{raw}': {e}"))?,
                        );
                    }
                    "--end" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--end",
                            "features tfbs-track-similarity",
                        )?;
                        state.end_arg = Some(
                            raw.parse::<usize>()
                                .map_err(|e| format!("Invalid --end value '{raw}': {e}"))?,
                        );
                    }
                    "--ranking-metric" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--ranking-metric",
                            "features tfbs-track-similarity",
                        )?;
                        ranking_metric = match raw.trim() {
                            "raw_pearson" => TfbsTrackSimilarityRankingMetric::RawPearson,
                            "smoothed_pearson" => TfbsTrackSimilarityRankingMetric::SmoothedPearson,
                            "raw_spearman" => TfbsTrackSimilarityRankingMetric::RawSpearman,
                            "smoothed_spearman" => {
                                TfbsTrackSimilarityRankingMetric::SmoothedSpearman
                            }
                            other => {
                                return Err(format!(
                                    "Unsupported --ranking-metric value '{other}' (expected raw_pearson, smoothed_pearson, raw_spearman, or smoothed_spearman)"
                                ));
                            }
                        };
                    }
                    "--score-kind" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--score-kind",
                            "features tfbs-track-similarity",
                        )?;
                        score_kind = match raw.trim() {
                            "llr_bits" => TfbsScoreTrackValueKind::LlrBits,
                            "llr_quantile" => TfbsScoreTrackValueKind::LlrQuantile,
                            "llr_background_quantile" => {
                                TfbsScoreTrackValueKind::LlrBackgroundQuantile
                            }
                            "llr_background_tail_log10" => {
                                TfbsScoreTrackValueKind::LlrBackgroundTailLog10
                            }
                            "true_log_odds_bits" => TfbsScoreTrackValueKind::TrueLogOddsBits,
                            "true_log_odds_quantile" => {
                                TfbsScoreTrackValueKind::TrueLogOddsQuantile
                            }
                            "true_log_odds_background_quantile" => {
                                TfbsScoreTrackValueKind::TrueLogOddsBackgroundQuantile
                            }
                            "true_log_odds_background_tail_log10" => {
                                TfbsScoreTrackValueKind::TrueLogOddsBackgroundTailLog10
                            }
                            other => {
                                return Err(format!(
                                    "Unsupported --score-kind value '{other}' (expected llr_bits, llr_quantile, llr_background_quantile, llr_background_tail_log10, true_log_odds_bits, true_log_odds_quantile, true_log_odds_background_quantile, or true_log_odds_background_tail_log10)"
                                ));
                            }
                        };
                    }
                    "--allow-negative" => {
                        clip_negative = false;
                        idx += 1;
                    }
                    "--species" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--species",
                            "features tfbs-track-similarity",
                        )?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            return Err("--species must not be empty".to_string());
                        }
                        species_filters.push(trimmed.to_string());
                    }
                    "--species-filter" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--species-filter",
                            "features tfbs-track-similarity",
                        )?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            return Err("--species must not be empty".to_string());
                        }
                        species_filters.push(trimmed.to_string());
                    }
                    "--include-remote-metadata" => {
                        include_remote_metadata = true;
                        idx += 1;
                    }
                    "--limit" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--limit",
                            "features tfbs-track-similarity",
                        )?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --limit value '{raw}': {e}"))?;
                        if parsed == 0 {
                            return Err("--limit must be >= 1".to_string());
                        }
                        limit = Some(parsed);
                    }
                    "--path" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--path",
                            "features tfbs-track-similarity",
                        )?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            path = None;
                        } else {
                            path = Some(trimmed.to_string());
                        }
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for features tfbs-track-similarity"
                        ));
                    }
                }
            }
            let anchor_motif = anchor_motif.ok_or_else(|| {
                "features tfbs-track-similarity requires --anchor-motif TOKEN".to_string()
            })?;
            let target = build_sequence_scan_target_from_feature_state(
                seq_id,
                sequence_text,
                topology,
                id_hint,
                state,
                "features tfbs-track-similarity",
            )?;
            Ok(ShellCommand::FeaturesTfbsTrackSimilarity {
                target,
                anchor_motif,
                candidate_motifs,
                ranking_metric,
                score_kind,
                clip_negative,
                species_filters,
                include_remote_metadata,
                limit,
                path,
            })
        }
        "tfbs-score-track-correlation-svg" => {
            if tokens.len() < 4 {
                return Err(
                    "features tfbs-score-track-correlation-svg requires SEQ_ID OUTPUT.svg --motif TOKEN [--motif TOKEN ...] [--range START..END|--start N --end N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--correlation-metric pearson|spearman] [--allow-negative]"
                        .to_string(),
                );
            }
            let seq_id = tokens[2].trim().to_string();
            if seq_id.is_empty() {
                return Err(
                    "features tfbs-score-track-correlation-svg requires non-empty SEQ_ID"
                        .to_string(),
                );
            }
            let output = tokens[3].clone();
            let mut motifs: Vec<String> = vec![];
            let mut start_0based: Option<usize> = None;
            let mut end_0based_exclusive: Option<usize> = None;
            let mut score_kind = TfbsScoreTrackValueKind::LlrBits;
            let mut correlation_metric = TfbsScoreTrackCorrelationMetric::Pearson;
            let mut correlation_signal_source = TfbsScoreTrackCorrelationSignalSource::MaxStrands;
            let mut clip_negative = true;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--motif" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--motif",
                            "features tfbs-score-track-correlation-svg",
                        )?;
                        let normalized = raw.trim();
                        if normalized.is_empty() {
                            return Err("--motif must not be empty".to_string());
                        }
                        motifs.push(normalized.to_string());
                    }
                    "--motifs" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--motifs",
                            "features tfbs-score-track-correlation-svg",
                        )?;
                        for token in raw
                            .split(',')
                            .map(str::trim)
                            .filter(|value| !value.is_empty())
                        {
                            motifs.push(token.to_string());
                        }
                    }
                    "--range" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--range",
                            "features tfbs-score-track-correlation-svg",
                        )?;
                        let (start, end) =
                            parse_feature_range(&raw, "features tfbs-score-track-correlation-svg")?;
                        start_0based = Some(start);
                        end_0based_exclusive = Some(end);
                    }
                    "--start" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--start",
                            "features tfbs-score-track-correlation-svg",
                        )?;
                        start_0based = Some(
                            raw.parse::<usize>()
                                .map_err(|e| format!("Invalid --start value '{raw}': {e}"))?,
                        );
                    }
                    "--end" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--end",
                            "features tfbs-score-track-correlation-svg",
                        )?;
                        end_0based_exclusive = Some(
                            raw.parse::<usize>()
                                .map_err(|e| format!("Invalid --end value '{raw}': {e}"))?,
                        );
                    }
                    "--score-kind" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--score-kind",
                            "features tfbs-score-track-correlation-svg",
                        )?;
                        score_kind = match raw.trim() {
                            "llr_bits" => TfbsScoreTrackValueKind::LlrBits,
                            "llr_quantile" => TfbsScoreTrackValueKind::LlrQuantile,
                            "llr_background_quantile" => {
                                TfbsScoreTrackValueKind::LlrBackgroundQuantile
                            }
                            "llr_background_tail_log10" => {
                                TfbsScoreTrackValueKind::LlrBackgroundTailLog10
                            }
                            "true_log_odds_bits" => TfbsScoreTrackValueKind::TrueLogOddsBits,
                            "true_log_odds_quantile" => {
                                TfbsScoreTrackValueKind::TrueLogOddsQuantile
                            }
                            "true_log_odds_background_quantile" => {
                                TfbsScoreTrackValueKind::TrueLogOddsBackgroundQuantile
                            }
                            "true_log_odds_background_tail_log10" => {
                                TfbsScoreTrackValueKind::TrueLogOddsBackgroundTailLog10
                            }
                            other => {
                                return Err(format!(
                                    "Unsupported --score-kind value '{other}' (expected llr_bits, llr_quantile, llr_background_quantile, llr_background_tail_log10, true_log_odds_bits, true_log_odds_quantile, true_log_odds_background_quantile, or true_log_odds_background_tail_log10)"
                                ));
                            }
                        };
                    }
                    "--correlation-metric" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--correlation-metric",
                            "features tfbs-score-track-correlation-svg",
                        )?;
                        correlation_metric = match raw.trim() {
                            "pearson" => TfbsScoreTrackCorrelationMetric::Pearson,
                            "spearman" => TfbsScoreTrackCorrelationMetric::Spearman,
                            other => {
                                return Err(format!(
                                    "Unsupported --correlation-metric value '{other}' (expected pearson or spearman)"
                                ));
                            }
                        };
                    }
                    "--signal-source" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--signal-source",
                            "features tfbs-score-track-correlation-svg",
                        )?;
                        correlation_signal_source = match raw.trim() {
                            "max_strands" => TfbsScoreTrackCorrelationSignalSource::MaxStrands,
                            "forward_only" => TfbsScoreTrackCorrelationSignalSource::ForwardOnly,
                            "reverse_only" => TfbsScoreTrackCorrelationSignalSource::ReverseOnly,
                            other => {
                                return Err(format!(
                                    "Unsupported --signal-source value '{other}' (expected max_strands, forward_only, or reverse_only)"
                                ));
                            }
                        };
                    }
                    "--allow-negative" => {
                        clip_negative = false;
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for features tfbs-score-track-correlation-svg"
                        ));
                    }
                }
            }
            if motifs.is_empty() {
                return Err(
                    "features tfbs-score-track-correlation-svg requires at least one --motif TOKEN"
                        .to_string(),
                );
            }
            let start_0based = start_0based.unwrap_or(0);
            let end_0based_exclusive = end_0based_exclusive.ok_or_else(|| {
                "features tfbs-score-track-correlation-svg requires --range START..END or --end N"
                    .to_string()
            })?;
            if start_0based >= end_0based_exclusive {
                return Err(format!(
                    "features tfbs-score-track-correlation-svg requires start < end (got {}..{})",
                    start_0based, end_0based_exclusive
                ));
            }
            Ok(ShellCommand::FeaturesTfbsScoreTrackCorrelationSvg {
                seq_id,
                motifs,
                start_0based,
                end_0based_exclusive,
                score_kind,
                correlation_metric,
                correlation_signal_source,
                clip_negative,
                output,
            })
        }
        "tfbs-scan" => {
            if tokens.len() < 3 {
                return Err(
                    "features tfbs-scan requires either SEQ_ID or --sequence-text DNA [--topology linear|circular] [--id-hint TEXT] --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--range START..END|--start N --end N] [--min-llr-bits VALUE] [--min-llr-quantile VALUE] [--per-tf-min-llr-bits TF=VALUE] [--per-tf-min-llr-quantile TF=VALUE] [--max-hits N] [--path FILE.json]"
                        .to_string(),
                );
            }
            let mut seq_id: Option<String> = None;
            let mut sequence_text: Option<String> = None;
            let mut topology = InlineSequenceTopology::Linear;
            let mut id_hint: Option<String> = None;
            let mut motifs: Vec<String> = vec![];
            let mut min_llr_bits: Option<f64> = None;
            let mut min_llr_quantile: Option<f64> = None;
            let mut per_tf_thresholds: std::collections::BTreeMap<String, TfThresholdOverride> =
                std::collections::BTreeMap::new();
            let mut max_hits: Option<usize> = None;
            let mut path: Option<String> = None;
            let mut state = FeatureQueryOptionState::default();
            let mut idx = 2usize;
            if tokens[2].starts_with("--") {
                // options-only inline form
            } else {
                let raw = tokens[2].trim().to_string();
                if !raw.is_empty() {
                    seq_id = Some(raw);
                }
                idx = 3usize;
            }
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--sequence-text" => {
                        sequence_text = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--sequence-text",
                            "features tfbs-scan",
                        )?);
                    }
                    "--topology" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--topology",
                            "features tfbs-scan",
                        )?;
                        topology = parse_inline_sequence_topology(&raw)?;
                    }
                    "--id-hint" => {
                        id_hint = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--id-hint",
                            "features tfbs-scan",
                        )?);
                    }
                    "--motif" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--motif", "features tfbs-scan")?;
                        let normalized = raw.trim();
                        if normalized.is_empty() {
                            return Err("--motif must not be empty".to_string());
                        }
                        motifs.push(normalized.to_string());
                    }
                    "--motifs" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--motifs", "features tfbs-scan")?;
                        for token in raw
                            .split(',')
                            .map(str::trim)
                            .filter(|value| !value.is_empty())
                        {
                            motifs.push(token.to_string());
                        }
                    }
                    "--min-llr-bits" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-llr-bits",
                            "features tfbs-scan",
                        )?;
                        min_llr_bits = Some(raw.parse::<f64>().map_err(|e| {
                            format!(
                                "Invalid --min-llr-bits value '{raw}' for features tfbs-scan: {e}"
                            )
                        })?);
                    }
                    "--min-llr-quantile" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-llr-quantile",
                            "features tfbs-scan",
                        )?;
                        min_llr_quantile = Some(raw.parse::<f64>().map_err(|e| {
                            format!(
                                "Invalid --min-llr-quantile value '{raw}' for features tfbs-scan: {e}"
                            )
                        })?);
                    }
                    "--per-tf-min-llr-bits" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--per-tf-min-llr-bits",
                            "features tfbs-scan",
                        )?;
                        let (tf, value) = parse_tf_threshold_override_value(
                            &raw,
                            "--per-tf-min-llr-bits",
                            "features tfbs-scan",
                        )?;
                        per_tf_thresholds
                            .entry(tf.clone())
                            .or_insert_with(|| TfThresholdOverride {
                                tf: tf.clone(),
                                min_llr_bits: None,
                                min_llr_quantile: None,
                            })
                            .min_llr_bits = Some(value);
                    }
                    "--per-tf-min-llr-quantile" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--per-tf-min-llr-quantile",
                            "features tfbs-scan",
                        )?;
                        let (tf, value) = parse_tf_threshold_override_value(
                            &raw,
                            "--per-tf-min-llr-quantile",
                            "features tfbs-scan",
                        )?;
                        per_tf_thresholds
                            .entry(tf.clone())
                            .or_insert_with(|| TfThresholdOverride {
                                tf: tf.clone(),
                                min_llr_bits: None,
                                min_llr_quantile: None,
                            })
                            .min_llr_quantile = Some(value);
                    }
                    "--max-hits" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-hits",
                            "features tfbs-scan",
                        )?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --max-hits value '{raw}' for features tfbs-scan: {e}")
                        })?;
                        if parsed == 0 {
                            return Err("--max-hits must be >= 1".to_string());
                        }
                        max_hits = Some(parsed);
                    }
                    "--path" => {
                        path = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--path",
                            "features tfbs-scan",
                        )?);
                    }
                    "--range" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--range", "features tfbs-scan")?;
                        if state.range_arg.is_some() {
                            return Err("--range was specified multiple times".to_string());
                        }
                        state.range_arg = Some(parse_feature_range(&raw, "features tfbs-scan")?);
                    }
                    "--start" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--start", "features tfbs-scan")?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --start value '{raw}' for features tfbs-scan: {e}")
                        })?;
                        state.start_arg = Some(parsed);
                    }
                    "--end" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--end", "features tfbs-scan")?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --end value '{raw}' for features tfbs-scan: {e}")
                        })?;
                        state.end_arg = Some(parsed);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for features tfbs-scan"));
                    }
                }
            }
            if motifs.is_empty() {
                return Err("features tfbs-scan requires at least one --motif TOKEN".to_string());
            }
            let target = build_sequence_scan_target_from_feature_state(
                seq_id,
                sequence_text,
                topology,
                id_hint,
                state,
                "features tfbs-scan",
            )?;
            Ok(ShellCommand::FeaturesTfbsScan {
                target,
                motifs,
                min_llr_bits,
                min_llr_quantile,
                per_tf_thresholds: per_tf_thresholds.into_values().collect(),
                max_hits,
                path,
            })
        }
        "restriction-scan" => {
            if tokens.len() < 3 {
                return Err(
                    "features restriction-scan requires either SEQ_ID or --sequence-text DNA [--topology linear|circular] [--id-hint TEXT] [--range START..END|--start N --end N] [--enzyme NAME] [--max-sites-per-enzyme N] [--no-cut-geometry] [--path FILE.json]"
                        .to_string(),
                );
            }
            let mut seq_id: Option<String> = None;
            let mut sequence_text: Option<String> = None;
            let mut topology = InlineSequenceTopology::Linear;
            let mut id_hint: Option<String> = None;
            let mut enzymes: Vec<String> = vec![];
            let mut max_sites_per_enzyme: Option<usize> = None;
            let mut include_cut_geometry = true;
            let mut path: Option<String> = None;
            let mut state = FeatureQueryOptionState::default();
            let mut idx = 2usize;
            if tokens[2].starts_with("--") {
                // options-only inline form
            } else {
                let raw = tokens[2].trim().to_string();
                if !raw.is_empty() {
                    seq_id = Some(raw);
                }
                idx = 3usize;
            }
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--sequence-text" => {
                        sequence_text = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--sequence-text",
                            "features restriction-scan",
                        )?);
                    }
                    "--topology" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--topology",
                            "features restriction-scan",
                        )?;
                        topology = parse_inline_sequence_topology(&raw)?;
                    }
                    "--id-hint" => {
                        id_hint = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--id-hint",
                            "features restriction-scan",
                        )?);
                    }
                    "--enzyme" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--enzyme",
                            "features restriction-scan",
                        )?;
                        let trimmed = raw.trim();
                        if !trimmed.is_empty() {
                            enzymes.push(trimmed.to_string());
                        }
                    }
                    "--max-sites-per-enzyme" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-sites-per-enzyme",
                            "features restriction-scan",
                        )?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --max-sites-per-enzyme value '{raw}' for features restriction-scan: {e}"
                            )
                        })?;
                        if parsed == 0 {
                            return Err("--max-sites-per-enzyme must be >= 1".to_string());
                        }
                        max_sites_per_enzyme = Some(parsed);
                    }
                    "--no-cut-geometry" => {
                        include_cut_geometry = false;
                        idx += 1;
                    }
                    "--path" => {
                        path = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--path",
                            "features restriction-scan",
                        )?);
                    }
                    "--range" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--range",
                            "features restriction-scan",
                        )?;
                        if state.range_arg.is_some() {
                            return Err("--range was specified multiple times".to_string());
                        }
                        state.range_arg =
                            Some(parse_feature_range(&raw, "features restriction-scan")?);
                    }
                    "--start" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--start",
                            "features restriction-scan",
                        )?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --start value '{raw}' for features restriction-scan: {e}"
                            )
                        })?;
                        state.start_arg = Some(parsed);
                    }
                    "--end" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--end",
                            "features restriction-scan",
                        )?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --end value '{raw}' for features restriction-scan: {e}"
                            )
                        })?;
                        state.end_arg = Some(parsed);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for features restriction-scan"
                        ));
                    }
                }
            }
            let target = build_sequence_scan_target_from_feature_state(
                seq_id,
                sequence_text,
                topology,
                id_hint,
                state,
                "features restriction-scan",
            )?;
            Ok(ShellCommand::FeaturesRestrictionScan {
                target,
                enzymes,
                max_sites_per_enzyme,
                include_cut_geometry,
                path,
            })
        }
        other => Err(format!(
            "Unknown features subcommand '{other}' (expected query, export-bed, tfbs-summary, tfbs-score-tracks-svg, tfbs-score-track-correlation-svg, tfbs-scan, or restriction-scan)"
        )),
    }
}

pub(super) fn build_seeded_primer_pair_operation(
    template: &str,
    roi_start_0based: usize,
    roi_end_0based_exclusive: usize,
) -> Operation {
    Operation::DesignPrimerPairs {
        template: template.to_string(),
        roi_start_0based,
        roi_end_0based: roi_end_0based_exclusive,
        forward: PrimerDesignSideConstraint::default(),
        reverse: PrimerDesignSideConstraint::default(),
        min_amplicon_bp: 120,
        max_amplicon_bp: 1200,
        pair_constraints: PrimerDesignPairConstraint::default(),
        max_tm_delta_c: Some(2.0),
        max_pairs: Some(200),
        report_id: None,
    }
}

pub(super) fn build_seeded_qpcr_operation(
    template: &str,
    roi_start_0based: usize,
    roi_end_0based_exclusive: usize,
) -> Operation {
    Operation::DesignQpcrAssays {
        template: template.to_string(),
        roi_start_0based,
        roi_end_0based: roi_end_0based_exclusive,
        forward: PrimerDesignSideConstraint::default(),
        reverse: PrimerDesignSideConstraint::default(),
        probe: PrimerDesignSideConstraint::default(),
        min_amplicon_bp: 120,
        max_amplicon_bp: 1200,
        pair_constraints: PrimerDesignPairConstraint::default(),
        max_tm_delta_c: Some(2.0),
        max_probe_tm_delta_c: Some(10.0),
        max_assays: Some(200),
        report_id: None,
    }
}

fn parse_restriction_cloning_handoff_mode(
    raw: &str,
) -> Result<RestrictionCloningPcrHandoffMode, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "single_site" | "single-site" | "single" => {
            Ok(RestrictionCloningPcrHandoffMode::SingleSite)
        }
        "directed_pair" | "directed-pair" | "directed" | "pair" => {
            Ok(RestrictionCloningPcrHandoffMode::DirectedPair)
        }
        other => Err(format!(
            "Unsupported restriction-cloning handoff mode '{other}', expected single_site|directed_pair"
        )),
    }
}

pub(super) fn parse_primers_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "primers requires a subcommand: design, design-qpcr, prepare-restriction-cloning, seed-restriction-cloning-handoff, restriction-cloning-vector-suggestions, list-restriction-cloning-handoffs, show-restriction-cloning-handoff, export-restriction-cloning-handoff, preflight, seed-from-feature, seed-from-splicing, seed-qpcr-from-feature, seed-qpcr-from-splicing, list-reports, show-report, export-report, list-qpcr-reports, show-qpcr-report, export-qpcr-report"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "design" => {
            if tokens.len() < 3 {
                return Err(
                    "primers design requires REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]"
                        .to_string(),
                );
            }
            let request_json = tokens[2].clone();
            let mut backend: Option<PrimerDesignBackend> = None;
            let mut primer3_executable: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--backend" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--backend", "primers design")?;
                        backend = Some(parse_primer_design_backend(&raw)?);
                    }
                    "--primer3-exec" | "--primer3-executable" => {
                        let flag = tokens[idx].clone();
                        primer3_executable = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            &flag,
                            "primers design",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for primers design"));
                    }
                }
            }
            Ok(ShellCommand::PrimersDesign {
                request_json,
                backend,
                primer3_executable,
            })
        }
        "design-qpcr" => {
            if tokens.len() < 3 {
                return Err(
                    "primers design-qpcr requires REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]"
                        .to_string(),
                );
            }
            let request_json = tokens[2].clone();
            let mut backend: Option<PrimerDesignBackend> = None;
            let mut primer3_executable: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--backend" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--backend",
                            "primers design-qpcr",
                        )?;
                        backend = Some(parse_primer_design_backend(&raw)?);
                    }
                    "--primer3-exec" | "--primer3-executable" => {
                        let flag = tokens[idx].clone();
                        primer3_executable = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            &flag,
                            "primers design-qpcr",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for primers design-qpcr"));
                    }
                }
            }
            Ok(ShellCommand::PrimersDesignQpcr {
                request_json,
                backend,
                primer3_executable,
            })
        }
        "prepare-restriction-cloning" => {
            if tokens.len() != 3 {
                return Err(
                    "primers prepare-restriction-cloning requires REQUEST_JSON_OR_@FILE"
                        .to_string(),
                );
            }
            Ok(ShellCommand::PrimersPrepareRestrictionCloning {
                request_json: tokens[2].clone(),
            })
        }
        "seed-restriction-cloning-handoff" => {
            if tokens.len() < 4 {
                return Err(
                    "primers seed-restriction-cloning-handoff requires PRIMER_REPORT_ID VECTOR_SEQ_ID [--pair-rank N] [--mode single_site|directed_pair] [--forward-enzyme NAME] [--reverse-enzyme NAME] [--forward-leader SEQ] [--reverse-leader SEQ]"
                        .to_string(),
                );
            }
            let primer_report_id = tokens[2].clone();
            let destination_vector_seq_id = tokens[3].clone();
            let mut pair_rank_1based = 1usize;
            let mut mode = RestrictionCloningPcrHandoffMode::SingleSite;
            let mut forward_enzyme: Option<String> = None;
            let mut reverse_enzyme: Option<String> = None;
            let mut forward_leader_5prime: Option<String> = None;
            let mut reverse_leader_5prime: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--pair-rank" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--pair-rank",
                            "primers seed-restriction-cloning-handoff",
                        )?;
                        pair_rank_1based = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --pair-rank value '{raw}': {e}"))?;
                        if pair_rank_1based == 0 {
                            return Err("--pair-rank must be >= 1".to_string());
                        }
                    }
                    "--mode" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--mode",
                            "primers seed-restriction-cloning-handoff",
                        )?;
                        mode = parse_restriction_cloning_handoff_mode(&raw)?;
                    }
                    "--forward-enzyme" => {
                        forward_enzyme = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--forward-enzyme",
                            "primers seed-restriction-cloning-handoff",
                        )?);
                    }
                    "--reverse-enzyme" => {
                        reverse_enzyme = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--reverse-enzyme",
                            "primers seed-restriction-cloning-handoff",
                        )?);
                    }
                    "--forward-leader" => {
                        forward_leader_5prime = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--forward-leader",
                            "primers seed-restriction-cloning-handoff",
                        )?);
                    }
                    "--reverse-leader" => {
                        reverse_leader_5prime = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--reverse-leader",
                            "primers seed-restriction-cloning-handoff",
                        )?);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for primers seed-restriction-cloning-handoff"
                        ));
                    }
                }
            }
            Ok(ShellCommand::PrimersSeedRestrictionCloningHandoff {
                primer_report_id,
                destination_vector_seq_id,
                pair_rank_1based,
                mode,
                forward_enzyme,
                reverse_enzyme,
                forward_leader_5prime,
                reverse_leader_5prime,
            })
        }
        "restriction-cloning-vector-suggestions" => {
            if tokens.len() != 3 {
                return Err(
                    "primers restriction-cloning-vector-suggestions requires SEQ_ID".to_string(),
                );
            }
            Ok(ShellCommand::PrimersRestrictionCloningVectorSuggestions {
                seq_id: tokens[2].clone(),
            })
        }
        "list-restriction-cloning-handoffs" => {
            if tokens.len() != 2 {
                return Err(
                    "primers list-restriction-cloning-handoffs takes no extra arguments"
                        .to_string(),
                );
            }
            Ok(ShellCommand::PrimersListRestrictionCloningHandoffs)
        }
        "show-restriction-cloning-handoff" => {
            if tokens.len() != 3 {
                return Err(
                    "primers show-restriction-cloning-handoff requires REPORT_ID".to_string(),
                );
            }
            Ok(ShellCommand::PrimersShowRestrictionCloningHandoff {
                report_id: tokens[2].clone(),
            })
        }
        "export-restriction-cloning-handoff" => {
            if tokens.len() != 4 {
                return Err(
                    "primers export-restriction-cloning-handoff requires REPORT_ID OUTPUT.json"
                        .to_string(),
                );
            }
            Ok(ShellCommand::PrimersExportRestrictionCloningHandoff {
                report_id: tokens[2].clone(),
                path: tokens[3].clone(),
            })
        }
        "preflight" => {
            let mut backend: Option<PrimerDesignBackend> = None;
            let mut primer3_executable: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--backend" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--backend", "primers preflight")?;
                        backend = Some(parse_primer_design_backend(&raw)?);
                    }
                    "--primer3-exec" | "--primer3-executable" => {
                        let flag = tokens[idx].clone();
                        primer3_executable = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            &flag,
                            "primers preflight",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for primers preflight"));
                    }
                }
            }
            Ok(ShellCommand::PrimersPreflight {
                backend,
                primer3_executable,
            })
        }
        "seed-from-feature" => {
            if tokens.len() != 4 {
                return Err("primers seed-from-feature requires SEQ_ID FEATURE_ID".to_string());
            }
            let seq_id = tokens[2].clone();
            let feature_id = tokens[3].parse::<usize>().map_err(|e| {
                format!(
                    "Invalid feature id '{}' for primers seed-from-feature: {e}",
                    tokens[3]
                )
            })?;
            Ok(ShellCommand::PrimersSeedFromFeature { seq_id, feature_id })
        }
        "seed-from-splicing" => {
            if tokens.len() != 4 {
                return Err("primers seed-from-splicing requires SEQ_ID FEATURE_ID".to_string());
            }
            let seq_id = tokens[2].clone();
            let feature_id = tokens[3].parse::<usize>().map_err(|e| {
                format!(
                    "Invalid feature id '{}' for primers seed-from-splicing: {e}",
                    tokens[3]
                )
            })?;
            Ok(ShellCommand::PrimersSeedFromSplicing { seq_id, feature_id })
        }
        "seed-qpcr-from-feature" => {
            if tokens.len() != 4 {
                return Err("primers seed-qpcr-from-feature requires SEQ_ID FEATURE_ID".to_string());
            }
            let seq_id = tokens[2].clone();
            let feature_id = tokens[3].parse::<usize>().map_err(|e| {
                format!(
                    "Invalid feature id '{}' for primers seed-qpcr-from-feature: {e}",
                    tokens[3]
                )
            })?;
            Ok(ShellCommand::PrimersSeedQpcrFromFeature { seq_id, feature_id })
        }
        "seed-qpcr-from-splicing" => {
            if tokens.len() != 4 {
                return Err(
                    "primers seed-qpcr-from-splicing requires SEQ_ID FEATURE_ID".to_string()
                );
            }
            let seq_id = tokens[2].clone();
            let feature_id = tokens[3].parse::<usize>().map_err(|e| {
                format!(
                    "Invalid feature id '{}' for primers seed-qpcr-from-splicing: {e}",
                    tokens[3]
                )
            })?;
            Ok(ShellCommand::PrimersSeedQpcrFromSplicing { seq_id, feature_id })
        }
        "list-reports" => {
            if tokens.len() != 2 {
                return Err("primers list-reports takes no options".to_string());
            }
            Ok(ShellCommand::PrimersListReports)
        }
        "show-report" => {
            if tokens.len() != 3 {
                return Err("primers show-report requires REPORT_ID".to_string());
            }
            Ok(ShellCommand::PrimersShowReport {
                report_id: tokens[2].clone(),
            })
        }
        "export-report" => {
            if tokens.len() != 4 {
                return Err("primers export-report requires REPORT_ID OUTPUT.json".to_string());
            }
            Ok(ShellCommand::PrimersExportReport {
                report_id: tokens[2].clone(),
                path: tokens[3].clone(),
            })
        }
        "list-qpcr-reports" => {
            if tokens.len() != 2 {
                return Err("primers list-qpcr-reports takes no options".to_string());
            }
            Ok(ShellCommand::PrimersListQpcrReports)
        }
        "show-qpcr-report" => {
            if tokens.len() != 3 {
                return Err("primers show-qpcr-report requires REPORT_ID".to_string());
            }
            Ok(ShellCommand::PrimersShowQpcrReport {
                report_id: tokens[2].clone(),
            })
        }
        "export-qpcr-report" => {
            if tokens.len() != 4 {
                return Err("primers export-qpcr-report requires REPORT_ID OUTPUT.json".to_string());
            }
            Ok(ShellCommand::PrimersExportQpcrReport {
                report_id: tokens[2].clone(),
                path: tokens[3].clone(),
            })
        }
        other => Err(format!(
            "Unknown primers subcommand '{other}' (expected design, design-qpcr, prepare-restriction-cloning, seed-restriction-cloning-handoff, restriction-cloning-vector-suggestions, list-restriction-cloning-handoffs, show-restriction-cloning-handoff, export-restriction-cloning-handoff, preflight, seed-from-feature, seed-from-splicing, seed-qpcr-from-feature, seed-qpcr-from-splicing, list-reports, show-report, export-report, list-qpcr-reports, show-qpcr-report, export-qpcr-report)"
        )),
    }
}

pub(super) fn parse_dotplot_mode(raw: &str) -> Result<DotplotMode, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "self_forward" | "self-forward" | "forward" | "self" => Ok(DotplotMode::SelfForward),
        "self_reverse_complement"
        | "self-reverse-complement"
        | "self_revcomp"
        | "self-revcomp"
        | "revcomp"
        | "reverse_complement"
        | "reverse-complement" => Ok(DotplotMode::SelfReverseComplement),
        "pair_forward" | "pair-forward" | "pair" => Ok(DotplotMode::PairForward),
        "pair_reverse_complement" | "pair-reverse-complement" | "pair_revcomp" | "pair-revcomp" => {
            Ok(DotplotMode::PairReverseComplement)
        }
        other => Err(format!(
            "Unsupported dotplot mode '{other}', expected self_forward|self_reverse_complement|pair_forward|pair_reverse_complement"
        )),
    }
}

pub(super) fn parse_flexibility_model(raw: &str) -> Result<FlexibilityModel, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "at_richness" | "at-richness" | "at" | "at_content" | "at-content" => {
            Ok(FlexibilityModel::AtRichness)
        }
        "at_skew" | "at-skew" | "skew" => Ok(FlexibilityModel::AtSkew),
        other => Err(format!(
            "Unsupported flexibility model '{other}', expected at_richness or at_skew"
        )),
    }
}

pub(super) fn parse_pairwise_alignment_mode(raw: &str) -> Result<PairwiseAlignmentMode, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "global" => Ok(PairwiseAlignmentMode::Global),
        "local" => Ok(PairwiseAlignmentMode::Local),
        other => Err(format!(
            "Unsupported alignment mode '{other}', expected global or local"
        )),
    }
}

pub(super) fn parse_rna_read_profile(raw: &str) -> Result<RnaReadInterpretationProfile, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "nanopore_cdna_v1" | "nanopore" | "nanopore_cdna" => {
            Ok(RnaReadInterpretationProfile::NanoporeCdnaV1)
        }
        "short_read_v1" | "shortread" | "short_read" => {
            Ok(RnaReadInterpretationProfile::ShortReadV1)
        }
        "transposon_v1" | "transposon" => Ok(RnaReadInterpretationProfile::TransposonV1),
        other => Err(format!(
            "Unsupported RNA-read profile '{other}', expected nanopore_cdna_v1|short_read_v1|transposon_v1"
        )),
    }
}

pub(super) fn parse_rna_read_input_format(raw: &str) -> Result<RnaReadInputFormat, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "fasta" | "fa" => Ok(RnaReadInputFormat::Fasta),
        other => Err(format!(
            "Unsupported RNA-read input format '{other}', expected fasta"
        )),
    }
}

pub(super) fn parse_splicing_scope_preset(raw: &str) -> Result<SplicingScopePreset, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "all_overlapping_both_strands" | "all" | "broad" => {
            Ok(SplicingScopePreset::AllOverlappingBothStrands)
        }
        "target_group_any_strand" | "target_group" | "group" => {
            Ok(SplicingScopePreset::TargetGroupAnyStrand)
        }
        "all_overlapping_target_strand" | "target_strand" | "strand" => {
            Ok(SplicingScopePreset::AllOverlappingTargetStrand)
        }
        "target_group_target_strand" | "group_strand" | "legacy" => {
            Ok(SplicingScopePreset::TargetGroupTargetStrand)
        }
        other => Err(format!(
            "Unsupported scope preset '{other}', expected all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand"
        )),
    }
}

pub(super) fn parse_rna_read_origin_mode(raw: &str) -> Result<RnaReadOriginMode, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "single_gene" | "single-gene" | "single" | "legacy" => Ok(RnaReadOriginMode::SingleGene),
        "multi_gene_sparse" | "multi-gene-sparse" | "multi_sparse" | "multi" => {
            Ok(RnaReadOriginMode::MultiGeneSparse)
        }
        other => Err(format!(
            "Unsupported origin mode '{other}', expected single_gene|multi_gene_sparse"
        )),
    }
}

pub(super) fn parse_rna_read_report_mode(raw: &str) -> Result<RnaReadReportMode, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "full" => Ok(RnaReadReportMode::Full),
        "seed_passed_only" | "seed-passed-only" | "seed_passed" | "seed-only" => {
            Ok(RnaReadReportMode::SeedPassedOnly)
        }
        other => Err(format!(
            "Unsupported report mode '{other}', expected full|seed_passed_only"
        )),
    }
}

pub(super) fn parse_rna_read_hit_selection(raw: &str) -> Result<RnaReadHitSelection, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "all" => Ok(RnaReadHitSelection::All),
        "seed_passed" | "seed" => Ok(RnaReadHitSelection::SeedPassed),
        "aligned" => Ok(RnaReadHitSelection::Aligned),
        other => Err(format!(
            "Unsupported hit selection '{other}', expected all|seed_passed|aligned"
        )),
    }
}

pub(super) fn parse_rna_read_gene_support_complete_rule(
    raw: &str,
) -> Result<RnaReadGeneSupportCompleteRule, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "near" => Ok(RnaReadGeneSupportCompleteRule::Near),
        "strict" => Ok(RnaReadGeneSupportCompleteRule::Strict),
        "exact" => Ok(RnaReadGeneSupportCompleteRule::Exact),
        other => Err(format!(
            "Unsupported RNA-read gene-support complete rule '{other}', expected near|strict|exact"
        )),
    }
}

pub(super) fn parse_rna_read_gene_support_audit_cohort_filter(
    raw: &str,
) -> Result<RnaReadGeneSupportAuditCohortFilter, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "all" => Ok(RnaReadGeneSupportAuditCohortFilter::All),
        "accepted" => Ok(RnaReadGeneSupportAuditCohortFilter::Accepted),
        "fragment" => Ok(RnaReadGeneSupportAuditCohortFilter::Fragment),
        "complete" => Ok(RnaReadGeneSupportAuditCohortFilter::Complete),
        "rejected" => Ok(RnaReadGeneSupportAuditCohortFilter::Rejected),
        other => Err(format!(
            "Unsupported RNA-read gene-support cohort filter '{other}', expected all|accepted|fragment|complete|rejected"
        )),
    }
}

pub(super) fn parse_rna_read_alignment_effect_filter(
    raw: &str,
) -> Result<RnaReadAlignmentInspectionEffectFilter, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "all_aligned" | "all-aligned" | "all" => {
            Ok(RnaReadAlignmentInspectionEffectFilter::AllAligned)
        }
        "confirmed_only" | "confirmed-only" | "confirmed" => {
            Ok(RnaReadAlignmentInspectionEffectFilter::ConfirmedOnly)
        }
        "disagreement_only" | "disagreement-only" | "disagreement" => {
            Ok(RnaReadAlignmentInspectionEffectFilter::DisagreementOnly)
        }
        "reassigned_only" | "reassigned-only" | "reassigned" => {
            Ok(RnaReadAlignmentInspectionEffectFilter::ReassignedOnly)
        }
        "no_phase1_only" | "no-phase1-only" | "no_phase1" | "no-phase1" => {
            Ok(RnaReadAlignmentInspectionEffectFilter::NoPhase1Only)
        }
        "selected_only" | "selected-only" | "selected" => {
            Ok(RnaReadAlignmentInspectionEffectFilter::SelectedOnly)
        }
        other => Err(format!(
            "Unsupported RNA-read alignment effect filter '{other}', expected all_aligned|confirmed_only|disagreement_only|reassigned_only|no_phase1_only|selected_only"
        )),
    }
}

pub(super) fn parse_rna_read_alignment_sort_key(
    raw: &str,
) -> Result<RnaReadAlignmentInspectionSortKey, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "rank" => Ok(RnaReadAlignmentInspectionSortKey::Rank),
        "identity" => Ok(RnaReadAlignmentInspectionSortKey::Identity),
        "coverage" => Ok(RnaReadAlignmentInspectionSortKey::Coverage),
        "score" => Ok(RnaReadAlignmentInspectionSortKey::Score),
        other => Err(format!(
            "Unsupported RNA-read alignment sort key '{other}', expected rank|identity|coverage|score"
        )),
    }
}

pub(super) fn parse_rna_read_score_density_variant(
    raw: &str,
) -> Result<RnaReadScoreDensityVariant, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "all_scored" | "all-scored" | "all" => Ok(RnaReadScoreDensityVariant::AllScored),
        "composite_seed_gate"
        | "composite-seed-gate"
        | "composite_gate"
        | "composite-gate"
        | "composite" => Ok(RnaReadScoreDensityVariant::CompositeSeedGate),
        other => Err(format!(
            "Unsupported RNA-read score-density variant '{other}', expected all_scored|composite_seed_gate"
        )),
    }
}

pub(super) fn parse_rna_read_record_indices(raw: &str) -> Result<Vec<usize>, String> {
    let mut indices = raw
        .split(',')
        .map(str::trim)
        .filter(|token| !token.is_empty())
        .map(|token| {
            token.parse::<usize>().map_err(|e| {
                format!("Invalid record index '{token}' (expected comma-separated non-negative integers): {e}")
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    if indices.is_empty() {
        return Err(
            "Invalid --record-indices value: expected comma-separated non-negative integers"
                .to_string(),
        );
    }
    indices.sort_unstable();
    indices.dedup();
    Ok(indices)
}

pub(super) fn parse_rna_read_score_density_scale(
    raw: &str,
) -> Result<RnaReadScoreDensityScale, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "linear" | "lin" => Ok(RnaReadScoreDensityScale::Linear),
        "log" | "log1p" => Ok(RnaReadScoreDensityScale::Log),
        other => Err(format!(
            "Unsupported score-density scale '{other}', expected linear|log"
        )),
    }
}

pub(super) fn parse_transcripts_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("transcripts requires a subcommand: derive".to_string());
    }
    match tokens[1].as_str() {
        "derive" => {
            if tokens.len() < 3 {
                return Err(
                    "transcripts derive requires SEQ_ID [--feature-id N ...] [--scope SCOPE] [--output-prefix PREFIX]"
                        .to_string(),
                );
            }
            let seq_id = tokens[2].trim().to_string();
            if seq_id.is_empty() {
                return Err("transcripts derive SEQ_ID must not be empty".to_string());
            }
            let mut feature_ids: Vec<usize> = vec![];
            let mut scope: Option<SplicingScopePreset> = None;
            let mut output_prefix: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--feature-id" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-id",
                            "transcripts derive",
                        )?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --feature-id value '{}' for transcripts derive: {e}",
                                raw
                            )
                        })?;
                        feature_ids.push(parsed);
                    }
                    "--scope" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--scope", "transcripts derive")?;
                        scope = Some(parse_splicing_scope_preset(&raw)?);
                    }
                    "--output-prefix" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--output-prefix",
                            "transcripts derive",
                        )?;
                        let trimmed = raw.trim();
                        output_prefix = if trimmed.is_empty() {
                            None
                        } else {
                            Some(trimmed.to_string())
                        };
                    }
                    other => {
                        return Err(format!("Unknown option '{}' for transcripts derive", other));
                    }
                }
            }
            if scope.is_some() && feature_ids.len() != 1 {
                return Err(
                    "transcripts derive --scope requires exactly one --feature-id seed".to_string(),
                );
            }
            Ok(ShellCommand::TranscriptsDerive {
                seq_id,
                feature_ids,
                scope,
                output_prefix,
            })
        }
        other => Err(format!(
            "Unknown transcripts subcommand '{}' (expected derive)",
            other
        )),
    }
}

pub(super) fn parse_dotplot_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "dotplot requires a subcommand: compute, overlay-compute, list, show, render-svg"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "compute" => {
            if tokens.len() < 3 {
                return Err(
                    "dotplot compute requires SEQ_ID [--reference-seq REF_SEQ_ID] [--start N] [--end N] [--ref-start N] [--ref-end N] [--mode MODE] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]"
                        .to_string(),
                );
            }
            let seq_id = tokens[2].trim().to_string();
            if seq_id.is_empty() {
                return Err("dotplot compute SEQ_ID must not be empty".to_string());
            }
            let mut reference_seq_id: Option<String> = None;
            let mut span_start_0based: Option<usize> = None;
            let mut span_end_0based: Option<usize> = None;
            let mut reference_span_start_0based: Option<usize> = None;
            let mut reference_span_end_0based: Option<usize> = None;
            let mut mode = DotplotMode::SelfForward;
            let mut word_size = 12usize;
            let mut step_bp = 2usize;
            let mut max_mismatches = 0usize;
            let mut tile_bp: Option<usize> = None;
            let mut dotplot_id: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--reference-seq" | "--ref-seq" => {
                        let flag = tokens[idx].clone();
                        let raw = parse_option_path(tokens, &mut idx, &flag, "dotplot compute")?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            return Err(format!("{flag} requires a non-empty sequence id"));
                        }
                        reference_seq_id = Some(trimmed.to_string());
                    }
                    "--start" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--start", "dotplot compute")?;
                        span_start_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --start value '{raw}' for dotplot compute: {e}")
                        })?);
                    }
                    "--end" => {
                        let raw = parse_option_path(tokens, &mut idx, "--end", "dotplot compute")?;
                        span_end_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --end value '{raw}' for dotplot compute: {e}")
                        })?);
                    }
                    "--ref-start" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--ref-start", "dotplot compute")?;
                        reference_span_start_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --ref-start value '{raw}' for dotplot compute: {e}")
                        })?);
                    }
                    "--ref-end" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--ref-end", "dotplot compute")?;
                        reference_span_end_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --ref-end value '{raw}' for dotplot compute: {e}")
                        })?);
                    }
                    "--mode" => {
                        let raw = parse_option_path(tokens, &mut idx, "--mode", "dotplot compute")?;
                        mode = parse_dotplot_mode(&raw)?;
                    }
                    "--word-size" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--word-size", "dotplot compute")?;
                        word_size = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --word-size value '{raw}' for dotplot compute: {e}")
                        })?;
                    }
                    "--step" | "--step-bp" => {
                        let flag = tokens[idx].clone();
                        let raw = parse_option_path(tokens, &mut idx, &flag, "dotplot compute")?;
                        step_bp = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid {flag} value '{raw}' for dotplot compute: {e}")
                        })?;
                    }
                    "--max-mismatches" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-mismatches",
                            "dotplot compute",
                        )?;
                        max_mismatches = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --max-mismatches value '{raw}' for dotplot compute: {e}"
                            )
                        })?;
                    }
                    "--tile-bp" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--tile-bp", "dotplot compute")?;
                        tile_bp = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --tile-bp value '{raw}' for dotplot compute: {e}")
                        })?);
                    }
                    "--id" => {
                        let raw = parse_option_path(tokens, &mut idx, "--id", "dotplot compute")?;
                        dotplot_id = Some(raw);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for dotplot compute"));
                    }
                }
            }
            Ok(ShellCommand::DotplotCompute {
                seq_id,
                reference_seq_id,
                span_start_0based,
                span_end_0based,
                reference_span_start_0based,
                reference_span_end_0based,
                mode,
                word_size,
                step_bp,
                max_mismatches,
                tile_bp,
                dotplot_id,
            })
        }
        "overlay-compute" => {
            if tokens.len() < 3 {
                return Err(
                    "dotplot overlay-compute requires OWNER_SEQ_ID [--reference-seq REF_SEQ_ID] --query-spec JSON_OR_@FILE [--query-spec JSON_OR_@FILE ...] [--ref-start N] [--ref-end N] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]"
                        .to_string(),
                );
            }
            let owner_seq_id = tokens[2].trim().to_string();
            if owner_seq_id.is_empty() {
                return Err("dotplot overlay-compute OWNER_SEQ_ID must not be empty".to_string());
            }
            let mut reference_seq_id: Option<String> = None;
            let mut reference_span_start_0based: Option<usize> = None;
            let mut reference_span_end_0based: Option<usize> = None;
            let mut queries: Vec<DotplotOverlayQuerySpec> = vec![];
            let mut word_size = 12usize;
            let mut step_bp = 2usize;
            let mut max_mismatches = 0usize;
            let mut tile_bp: Option<usize> = None;
            let mut dotplot_id: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--reference-seq" | "--ref-seq" => {
                        let flag = tokens[idx].clone();
                        let raw =
                            parse_option_path(tokens, &mut idx, &flag, "dotplot overlay-compute")?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            return Err(format!("{flag} requires a non-empty sequence id"));
                        }
                        reference_seq_id = Some(trimmed.to_string());
                    }
                    "--query-spec" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--query-spec",
                            "dotplot overlay-compute",
                        )?;
                        let parsed = parse_required_json_payload::<DotplotOverlayQuerySpec>(
                            &raw,
                            "dotplot overlay query spec",
                        )?;
                        queries.push(parsed);
                    }
                    "--ref-start" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--ref-start",
                            "dotplot overlay-compute",
                        )?;
                        reference_span_start_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --ref-start value '{raw}' for dotplot overlay-compute: {e}"
                            )
                        })?);
                    }
                    "--ref-end" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--ref-end",
                            "dotplot overlay-compute",
                        )?;
                        reference_span_end_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --ref-end value '{raw}' for dotplot overlay-compute: {e}"
                            )
                        })?);
                    }
                    "--word-size" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--word-size",
                            "dotplot overlay-compute",
                        )?;
                        word_size = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --word-size value '{raw}' for dotplot overlay-compute: {e}"
                            )
                        })?;
                    }
                    "--step" | "--step-bp" => {
                        let flag = tokens[idx].clone();
                        let raw =
                            parse_option_path(tokens, &mut idx, &flag, "dotplot overlay-compute")?;
                        step_bp = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid {flag} value '{raw}' for dotplot overlay-compute: {e}")
                        })?;
                    }
                    "--max-mismatches" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-mismatches",
                            "dotplot overlay-compute",
                        )?;
                        max_mismatches = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --max-mismatches value '{raw}' for dotplot overlay-compute: {e}"
                            )
                        })?;
                    }
                    "--tile-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--tile-bp",
                            "dotplot overlay-compute",
                        )?;
                        tile_bp = Some(raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --tile-bp value '{raw}' for dotplot overlay-compute: {e}"
                            )
                        })?);
                    }
                    "--id" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--id", "dotplot overlay-compute")?;
                        dotplot_id = Some(raw);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for dotplot overlay-compute"
                        ));
                    }
                }
            }
            if queries.is_empty() {
                return Err(
                    "dotplot overlay-compute requires at least one --query-spec JSON_OR_@FILE"
                        .to_string(),
                );
            }
            Ok(ShellCommand::DotplotOverlayCompute {
                owner_seq_id,
                reference_seq_id,
                reference_span_start_0based,
                reference_span_end_0based,
                queries,
                word_size,
                step_bp,
                max_mismatches,
                tile_bp,
                dotplot_id,
            })
        }
        "list" => {
            if tokens.len() > 3 {
                return Err("dotplot list expects at most one optional SEQ_ID".to_string());
            }
            let seq_id = if tokens.len() == 3 {
                let value = tokens[2].trim();
                if value.is_empty() {
                    None
                } else {
                    Some(value.to_string())
                }
            } else {
                None
            };
            Ok(ShellCommand::DotplotList { seq_id })
        }
        "show" => {
            if tokens.len() != 3 {
                return Err("dotplot show requires DOTPLOT_ID".to_string());
            }
            Ok(ShellCommand::DotplotShow {
                dotplot_id: tokens[2].clone(),
            })
        }
        "render-svg" => {
            if tokens.len() < 5 {
                return Err(
                    "dotplot render-svg requires SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor|query_anchor_bp] [--overlay-anchor-exon START..END]"
                        .to_string(),
                );
            }
            let seq_id = tokens[2].trim();
            let dotplot_id = tokens[3].trim();
            if seq_id.is_empty() {
                return Err("dotplot render-svg requires non-empty SEQ_ID".to_string());
            }
            if dotplot_id.is_empty() {
                return Err("dotplot render-svg requires non-empty DOTPLOT_ID".to_string());
            }
            let output = tokens[4].clone();
            let mut flex_track_id: Option<String> = None;
            let mut display_density_threshold: Option<f32> = None;
            let mut display_intensity_gain: Option<f32> = None;
            let mut overlay_x_axis_mode = DotplotOverlayXAxisMode::PercentLength;
            let mut overlay_anchor_exon: Option<DotplotOverlayAnchorExonRef> = None;
            let mut idx = 5usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--flex-track" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--flex-track",
                            "dotplot render-svg",
                        )?;
                        let trimmed = raw.trim();
                        if !trimmed.is_empty() {
                            flex_track_id = Some(trimmed.to_string());
                        }
                    }
                    "--display-threshold" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--display-threshold",
                            "dotplot render-svg",
                        )?;
                        display_density_threshold = Some(raw.parse::<f32>().map_err(|e| {
                            format!(
                                "Invalid --display-threshold value '{raw}' for dotplot render-svg: {e}"
                            )
                        })?);
                    }
                    "--intensity-gain" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--intensity-gain",
                            "dotplot render-svg",
                        )?;
                        display_intensity_gain = Some(raw.parse::<f32>().map_err(|e| {
                            format!(
                                "Invalid --intensity-gain value '{raw}' for dotplot render-svg: {e}"
                            )
                        })?);
                    }
                    "--overlay-x-axis" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--overlay-x-axis",
                            "dotplot render-svg",
                        )?;
                        overlay_x_axis_mode = match raw.trim() {
                            "percent_length" => DotplotOverlayXAxisMode::PercentLength,
                            "left_aligned_bp" => DotplotOverlayXAxisMode::LeftAlignedBp,
                            "right_aligned_bp" => DotplotOverlayXAxisMode::RightAlignedBp,
                            "shared_exon_anchor" => DotplotOverlayXAxisMode::SharedExonAnchor,
                            "query_anchor_bp" => DotplotOverlayXAxisMode::QueryAnchorBp,
                            other => {
                                return Err(format!(
                                    "Invalid --overlay-x-axis value '{other}' for dotplot render-svg: expected percent_length, left_aligned_bp, right_aligned_bp, shared_exon_anchor, or query_anchor_bp"
                                ));
                            }
                        };
                    }
                    "--overlay-anchor-exon" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--overlay-anchor-exon",
                            "dotplot render-svg",
                        )?;
                        overlay_anchor_exon = Some(DotplotOverlayAnchorExonRef::parse(&raw)?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for dotplot render-svg"));
                    }
                }
            }
            Ok(ShellCommand::RenderDotplotSvg {
                seq_id: seq_id.to_string(),
                dotplot_id: dotplot_id.to_string(),
                output,
                flex_track_id,
                display_density_threshold,
                display_intensity_gain,
                overlay_x_axis_mode,
                overlay_anchor_exon,
            })
        }
        other => Err(format!(
            "Unknown dotplot subcommand '{other}' (expected compute, list, show, render-svg)"
        )),
    }
}

pub(super) fn parse_flex_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("flex requires a subcommand: compute, list, show".to_string());
    }
    match tokens[1].as_str() {
        "compute" => {
            if tokens.len() < 3 {
                return Err(
                    "flex compute requires SEQ_ID [--start N] [--end N] [--model MODEL] [--bin-bp N] [--smoothing-bp N] [--id TRACK_ID]"
                        .to_string(),
                );
            }
            let seq_id = tokens[2].trim().to_string();
            if seq_id.is_empty() {
                return Err("flex compute SEQ_ID must not be empty".to_string());
            }
            let mut span_start_0based: Option<usize> = None;
            let mut span_end_0based: Option<usize> = None;
            let mut model = FlexibilityModel::AtRichness;
            let mut bin_bp = 25usize;
            let mut smoothing_bp: Option<usize> = None;
            let mut track_id: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--start" => {
                        let raw = parse_option_path(tokens, &mut idx, "--start", "flex compute")?;
                        span_start_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --start value '{raw}' for flex compute: {e}")
                        })?);
                    }
                    "--end" => {
                        let raw = parse_option_path(tokens, &mut idx, "--end", "flex compute")?;
                        span_end_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --end value '{raw}' for flex compute: {e}")
                        })?);
                    }
                    "--model" => {
                        let raw = parse_option_path(tokens, &mut idx, "--model", "flex compute")?;
                        model = parse_flexibility_model(&raw)?;
                    }
                    "--bin-bp" => {
                        let raw = parse_option_path(tokens, &mut idx, "--bin-bp", "flex compute")?;
                        bin_bp = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --bin-bp value '{raw}' for flex compute: {e}")
                        })?;
                    }
                    "--smoothing-bp" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--smoothing-bp", "flex compute")?;
                        smoothing_bp = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --smoothing-bp value '{raw}' for flex compute: {e}")
                        })?);
                    }
                    "--id" => {
                        let raw = parse_option_path(tokens, &mut idx, "--id", "flex compute")?;
                        track_id = Some(raw);
                    }
                    other => return Err(format!("Unknown option '{other}' for flex compute")),
                }
            }
            Ok(ShellCommand::FlexCompute {
                seq_id,
                span_start_0based,
                span_end_0based,
                model,
                bin_bp,
                smoothing_bp,
                track_id,
            })
        }
        "list" => {
            if tokens.len() > 3 {
                return Err("flex list expects at most one optional SEQ_ID".to_string());
            }
            let seq_id = if tokens.len() == 3 {
                let value = tokens[2].trim();
                if value.is_empty() {
                    None
                } else {
                    Some(value.to_string())
                }
            } else {
                None
            };
            Ok(ShellCommand::FlexList { seq_id })
        }
        "show" => {
            if tokens.len() != 3 {
                return Err("flex show requires TRACK_ID".to_string());
            }
            Ok(ShellCommand::FlexShow {
                track_id: tokens[2].clone(),
            })
        }
        other => Err(format!(
            "Unknown flex subcommand '{other}' (expected compute, list, show)"
        )),
    }
}

pub(super) fn parse_splicing_refs_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("splicing-refs requires a subcommand: derive".to_string());
    }
    match tokens[1].as_str() {
        "derive" => {
            if tokens.len() < 5 {
                return Err(
                    "splicing-refs derive requires SEQ_ID START_0BASED END_0BASED [--seed-feature-id N] [--scope all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--output-prefix PREFIX]"
                        .to_string(),
                );
            }
            let seq_id = tokens[2].trim().to_string();
            if seq_id.is_empty() {
                return Err("splicing-refs derive SEQ_ID must not be empty".to_string());
            }
            let span_start_0based = tokens[3].parse::<usize>().map_err(|e| {
                format!(
                    "Invalid START_0BASED '{}' for splicing-refs derive: {e}",
                    tokens[3]
                )
            })?;
            let span_end_0based = tokens[4].parse::<usize>().map_err(|e| {
                format!(
                    "Invalid END_0BASED '{}' for splicing-refs derive: {e}",
                    tokens[4]
                )
            })?;
            let mut seed_feature_id: Option<usize> = None;
            let mut scope = SplicingScopePreset::TargetGroupTargetStrand;
            let mut output_prefix: Option<String> = None;
            let mut idx = 5usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--seed-feature-id" | "--seed-feature" => {
                        let flag = tokens[idx].clone();
                        let raw =
                            parse_option_path(tokens, &mut idx, &flag, "splicing-refs derive")?;
                        seed_feature_id = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid {flag} value '{raw}' for splicing-refs derive: {e}")
                        })?);
                    }
                    "--scope" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--scope", "splicing-refs derive")?;
                        scope = parse_splicing_scope_preset(&raw)?;
                    }
                    "--output-prefix" | "--prefix" => {
                        let flag = tokens[idx].clone();
                        let raw =
                            parse_option_path(tokens, &mut idx, &flag, "splicing-refs derive")?;
                        let trimmed = raw.trim();
                        output_prefix = if trimmed.is_empty() {
                            None
                        } else {
                            Some(trimmed.to_string())
                        };
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for splicing-refs derive"));
                    }
                }
            }
            Ok(ShellCommand::SplicingRefsDerive {
                seq_id,
                span_start_0based,
                span_end_0based,
                seed_feature_id,
                scope,
                output_prefix,
            })
        }
        other => Err(format!(
            "Unknown splicing-refs subcommand '{other}' (expected derive)"
        )),
    }
}

pub(super) fn parse_align_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("align requires a subcommand: compute".to_string());
    }
    match tokens[1].as_str() {
        "compute" => {
            if tokens.len() < 4 {
                return Err(
                    "align compute requires QUERY_SEQ_ID TARGET_SEQ_ID [--query-start N] [--query-end N] [--target-start N] [--target-end N] [--mode global|local] [--match N] [--mismatch N] [--gap-open N] [--gap-extend N]"
                        .to_string(),
                );
            }
            let query_seq_id = tokens[2].trim().to_string();
            if query_seq_id.is_empty() {
                return Err("align compute QUERY_SEQ_ID must not be empty".to_string());
            }
            let target_seq_id = tokens[3].trim().to_string();
            if target_seq_id.is_empty() {
                return Err("align compute TARGET_SEQ_ID must not be empty".to_string());
            }
            let mut query_span_start_0based: Option<usize> = None;
            let mut query_span_end_0based: Option<usize> = None;
            let mut target_span_start_0based: Option<usize> = None;
            let mut target_span_end_0based: Option<usize> = None;
            let mut mode = PairwiseAlignmentMode::Global;
            let mut match_score = 2i32;
            let mut mismatch_score = -3i32;
            let mut gap_open = -5i32;
            let mut gap_extend = -1i32;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--query-start" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--query-start", "align compute")?;
                        query_span_start_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --query-start value '{raw}' for align compute: {e}")
                        })?);
                    }
                    "--query-end" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--query-end", "align compute")?;
                        query_span_end_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --query-end value '{raw}' for align compute: {e}")
                        })?);
                    }
                    "--target-start" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--target-start", "align compute")?;
                        target_span_start_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --target-start value '{raw}' for align compute: {e}")
                        })?);
                    }
                    "--target-end" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--target-end", "align compute")?;
                        target_span_end_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --target-end value '{raw}' for align compute: {e}")
                        })?);
                    }
                    "--mode" => {
                        let raw = parse_option_path(tokens, &mut idx, "--mode", "align compute")?;
                        mode = parse_pairwise_alignment_mode(&raw)?;
                    }
                    "--match" | "--match-score" => {
                        let flag = tokens[idx].clone();
                        let raw = parse_option_path(tokens, &mut idx, &flag, "align compute")?;
                        match_score = raw.parse::<i32>().map_err(|e| {
                            format!("Invalid {flag} value '{raw}' for align compute: {e}")
                        })?;
                    }
                    "--mismatch" | "--mismatch-score" => {
                        let flag = tokens[idx].clone();
                        let raw = parse_option_path(tokens, &mut idx, &flag, "align compute")?;
                        mismatch_score = raw.parse::<i32>().map_err(|e| {
                            format!("Invalid {flag} value '{raw}' for align compute: {e}")
                        })?;
                    }
                    "--gap-open" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--gap-open", "align compute")?;
                        gap_open = raw.parse::<i32>().map_err(|e| {
                            format!("Invalid --gap-open value '{raw}' for align compute: {e}")
                        })?;
                    }
                    "--gap-extend" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--gap-extend", "align compute")?;
                        gap_extend = raw.parse::<i32>().map_err(|e| {
                            format!("Invalid --gap-extend value '{raw}' for align compute: {e}")
                        })?;
                    }
                    other => return Err(format!("Unknown option '{other}' for align compute")),
                }
            }
            Ok(ShellCommand::AlignCompute {
                query_seq_id,
                target_seq_id,
                query_span_start_0based,
                query_span_end_0based,
                target_span_start_0based,
                target_span_end_0based,
                mode,
                match_score,
                mismatch_score,
                gap_open,
                gap_extend,
            })
        }
        other => Err(format!(
            "Unknown align subcommand '{other}' (expected compute)"
        )),
    }
}

pub(super) fn parse_seq_trace_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("seq-trace requires a subcommand: import, list, show".to_string());
    }
    match tokens[1].as_str() {
        "import" => {
            if tokens.len() < 3 {
                return Err(
                    "seq-trace import requires PATH [--trace-id ID] [--seq-id ID]".to_string(),
                );
            }
            let path = tokens[2].clone();
            let mut trace_id: Option<String> = None;
            let mut seq_id: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--trace-id" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--trace-id", "seq-trace import")?;
                        let trimmed = raw.trim();
                        if !trimmed.is_empty() {
                            trace_id = Some(trimmed.to_string());
                        }
                    }
                    "--seq-id" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--seq-id", "seq-trace import")?;
                        let trimmed = raw.trim();
                        if !trimmed.is_empty() {
                            seq_id = Some(trimmed.to_string());
                        }
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for seq-trace import"));
                    }
                }
            }
            Ok(ShellCommand::SeqTraceImport {
                path,
                trace_id,
                seq_id,
            })
        }
        "list" => {
            let mut seq_id: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--seq-id" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--seq-id", "seq-trace list")?;
                        let trimmed = raw.trim();
                        seq_id = if trimmed.is_empty() {
                            None
                        } else {
                            Some(trimmed.to_string())
                        };
                    }
                    other => {
                        if seq_id.is_none() && !other.starts_with('-') {
                            seq_id = Some(other.to_string());
                            idx += 1;
                        } else {
                            return Err(format!("Unknown option '{other}' for seq-trace list"));
                        }
                    }
                }
            }
            Ok(ShellCommand::SeqTraceList { seq_id })
        }
        "show" => {
            if tokens.len() != 3 {
                return Err("seq-trace show requires TRACE_ID".to_string());
            }
            Ok(ShellCommand::SeqTraceShow {
                trace_id: tokens[2].trim().to_string(),
            })
        }
        other => Err(format!(
            "Unknown seq-trace subcommand '{other}' (expected import, list, show)"
        )),
    }
}

pub(super) fn parse_seq_confirm_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "seq-confirm requires a subcommand: run, list-reports, show-report, export-report, export-support-tsv"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "run" => {
            if tokens.len() < 4 {
                return Err(
                    "seq-confirm run requires EXPECTED_SEQ_ID plus at least one --reads/--read or --trace-ids/--trace-id value [--baseline BASELINE_SEQ_ID] [--junction LEFT_END_0BASED]... [--junction-flank N] [--report-id ID] [--mode global|local] [--match N] [--mismatch N] [--gap-open N] [--gap-extend N] [--min-identity F] [--min-target-coverage F] [--allow-reverse-complement|--no-reverse-complement]"
                        .to_string(),
                );
            }
            let expected_seq_id = tokens[2].trim().to_string();
            if expected_seq_id.is_empty() {
                return Err("seq-confirm run EXPECTED_SEQ_ID must not be empty".to_string());
            }
            let mut baseline_seq_id: Option<String> = None;
            let mut read_seq_ids: Vec<String> = vec![];
            let mut trace_ids: Vec<String> = vec![];
            let mut targets: Vec<SequencingConfirmationTargetSpec> = vec![];
            let mut alignment_mode = PairwiseAlignmentMode::Local;
            let mut match_score = 2i32;
            let mut mismatch_score = -3i32;
            let mut gap_open = -5i32;
            let mut gap_extend = -1i32;
            let mut min_identity_fraction = 0.80f64;
            let mut min_target_coverage_fraction = 1.0f64;
            let mut allow_reverse_complement = true;
            let mut report_id: Option<String> = None;
            let mut junction_flank = 12usize;
            let mut junction_positions: Vec<usize> = vec![];
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--baseline" | "--baseline-seq-id" => {
                        let flag = tokens[idx].clone();
                        let raw = parse_option_path(tokens, &mut idx, &flag, "seq-confirm run")?;
                        let trimmed = raw.trim();
                        baseline_seq_id = if trimmed.is_empty() {
                            None
                        } else {
                            Some(trimmed.to_string())
                        };
                    }
                    "--reads" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--reads", "seq-confirm run")?;
                        let parsed = raw
                            .split(',')
                            .map(str::trim)
                            .filter(|value| !value.is_empty())
                            .map(|value| value.to_string())
                            .collect::<Vec<_>>();
                        if parsed.is_empty() {
                            return Err("--reads for seq-confirm run must include at least one ID"
                                .to_string());
                        }
                        read_seq_ids.extend(parsed);
                    }
                    "--read" => {
                        let raw = parse_option_path(tokens, &mut idx, "--read", "seq-confirm run")?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            return Err("--read for seq-confirm run must not be empty".to_string());
                        }
                        read_seq_ids.push(trimmed.to_string());
                    }
                    // Keep imported trace evidence flags symmetrical with the
                    // called-read flags so CLI shell and GUI shell stay in
                    // lockstep as trace-aware confirmation expands.
                    "--trace-ids" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--trace-ids", "seq-confirm run")?;
                        let parsed = raw
                            .split(',')
                            .map(str::trim)
                            .filter(|value| !value.is_empty())
                            .map(|value| value.to_string())
                            .collect::<Vec<_>>();
                        if parsed.is_empty() {
                            return Err(
                                "--trace-ids for seq-confirm run must include at least one ID"
                                    .to_string(),
                            );
                        }
                        trace_ids.extend(parsed);
                    }
                    "--trace-id" | "--trace" => {
                        let flag = tokens[idx].clone();
                        let raw = parse_option_path(tokens, &mut idx, &flag, "seq-confirm run")?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            return Err(format!("{flag} for seq-confirm run must not be empty"));
                        }
                        trace_ids.push(trimmed.to_string());
                    }
                    "--junction" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--junction", "seq-confirm run")?;
                        let left_end = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --junction value '{raw}' for seq-confirm run: {e}")
                        })?;
                        junction_positions.push(left_end);
                    }
                    "--junction-flank" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--junction-flank",
                            "seq-confirm run",
                        )?;
                        junction_flank = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --junction-flank value '{raw}' for seq-confirm run: {e}"
                            )
                        })?;
                        if junction_flank == 0 {
                            return Err(
                                "--junction-flank for seq-confirm run must be >= 1".to_string()
                            );
                        }
                    }
                    "--report-id" => {
                        report_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--report-id",
                            "seq-confirm run",
                        )?);
                    }
                    "--mode" => {
                        let raw = parse_option_path(tokens, &mut idx, "--mode", "seq-confirm run")?;
                        alignment_mode = parse_pairwise_alignment_mode(&raw)?;
                    }
                    "--match" | "--match-score" => {
                        let flag = tokens[idx].clone();
                        let raw = parse_option_path(tokens, &mut idx, &flag, "seq-confirm run")?;
                        match_score = raw.parse::<i32>().map_err(|e| {
                            format!("Invalid {flag} value '{raw}' for seq-confirm run: {e}")
                        })?;
                    }
                    "--mismatch" | "--mismatch-score" => {
                        let flag = tokens[idx].clone();
                        let raw = parse_option_path(tokens, &mut idx, &flag, "seq-confirm run")?;
                        mismatch_score = raw.parse::<i32>().map_err(|e| {
                            format!("Invalid {flag} value '{raw}' for seq-confirm run: {e}")
                        })?;
                    }
                    "--gap-open" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--gap-open", "seq-confirm run")?;
                        gap_open = raw.parse::<i32>().map_err(|e| {
                            format!("Invalid --gap-open value '{raw}' for seq-confirm run: {e}")
                        })?;
                    }
                    "--gap-extend" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--gap-extend", "seq-confirm run")?;
                        gap_extend = raw.parse::<i32>().map_err(|e| {
                            format!("Invalid --gap-extend value '{raw}' for seq-confirm run: {e}")
                        })?;
                    }
                    "--min-identity" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-identity",
                            "seq-confirm run",
                        )?;
                        min_identity_fraction = raw.parse::<f64>().map_err(|e| {
                            format!("Invalid --min-identity value '{raw}' for seq-confirm run: {e}")
                        })?;
                    }
                    "--min-target-coverage" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-target-coverage",
                            "seq-confirm run",
                        )?;
                        min_target_coverage_fraction = raw.parse::<f64>().map_err(|e| {
                            format!(
                                "Invalid --min-target-coverage value '{raw}' for seq-confirm run: {e}"
                            )
                        })?;
                    }
                    "--allow-reverse-complement" => {
                        allow_reverse_complement = true;
                        idx += 1;
                    }
                    "--no-reverse-complement" => {
                        allow_reverse_complement = false;
                        idx += 1;
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for seq-confirm run"));
                    }
                }
            }
            if read_seq_ids.is_empty() && trace_ids.is_empty() {
                return Err(
                    "seq-confirm run requires at least one read via --reads/--read or one trace via --trace-ids/--trace-id"
                        .to_string(),
                );
            }
            if !junction_positions.is_empty() {
                targets.extend(
                    junction_positions
                        .iter()
                        .enumerate()
                        .map(|(idx, left_end)| SequencingConfirmationTargetSpec {
                            target_id: format!("junction_{}", idx + 1),
                            label: format!("Junction @ {left_end}"),
                            kind: SequencingConfirmationTargetKind::Junction,
                            start_0based: left_end.saturating_sub(junction_flank),
                            end_0based_exclusive: left_end.saturating_add(junction_flank),
                            junction_left_end_0based: Some(*left_end),
                            expected_bases: None,
                            baseline_bases: None,
                            required: true,
                        }),
                );
            }
            Ok(ShellCommand::SeqConfirmRun {
                expected_seq_id,
                baseline_seq_id,
                read_seq_ids,
                trace_ids,
                targets,
                alignment_mode,
                match_score,
                mismatch_score,
                gap_open,
                gap_extend,
                min_identity_fraction,
                min_target_coverage_fraction,
                allow_reverse_complement,
                report_id,
            })
        }
        "list-reports" => {
            let mut expected_seq_id: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--expected" | "--expected-seq-id" => {
                        let flag = tokens[idx].clone();
                        let raw =
                            parse_option_path(tokens, &mut idx, &flag, "seq-confirm list-reports")?;
                        let trimmed = raw.trim();
                        expected_seq_id = if trimmed.is_empty() {
                            None
                        } else {
                            Some(trimmed.to_string())
                        };
                    }
                    other => {
                        if expected_seq_id.is_none() && !other.starts_with('-') {
                            expected_seq_id = Some(other.to_string());
                            idx += 1;
                        } else {
                            return Err(format!(
                                "Unknown option '{other}' for seq-confirm list-reports"
                            ));
                        }
                    }
                }
            }
            Ok(ShellCommand::SeqConfirmListReports { expected_seq_id })
        }
        "show-report" => {
            if tokens.len() != 3 {
                return Err("seq-confirm show-report requires REPORT_ID".to_string());
            }
            Ok(ShellCommand::SeqConfirmShowReport {
                report_id: tokens[2].trim().to_string(),
            })
        }
        "export-report" => {
            if tokens.len() != 4 {
                return Err("seq-confirm export-report requires REPORT_ID OUTPUT.json".to_string());
            }
            Ok(ShellCommand::SeqConfirmExportReport {
                report_id: tokens[2].trim().to_string(),
                path: tokens[3].clone(),
            })
        }
        "export-support-tsv" => {
            if tokens.len() != 4 {
                return Err(
                    "seq-confirm export-support-tsv requires REPORT_ID OUTPUT.tsv".to_string(),
                );
            }
            Ok(ShellCommand::SeqConfirmExportSupportTsv {
                report_id: tokens[2].trim().to_string(),
                path: tokens[3].clone(),
            })
        }
        other => Err(format!(
            "Unknown seq-confirm subcommand '{other}' (expected run, list-reports, show-report, export-report, export-support-tsv)"
        )),
    }
}

fn parse_translation_speed_profile_for_shell(
    raw: &str,
    context: &str,
) -> Result<TranslationSpeedProfile, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "human" => Ok(TranslationSpeedProfile::Human),
        "mouse" => Ok(TranslationSpeedProfile::Mouse),
        "yeast" => Ok(TranslationSpeedProfile::Yeast),
        "ecoli" | "e_coli" | "e-coli" => Ok(TranslationSpeedProfile::Ecoli),
        other => Err(format!(
            "Unknown --speed-profile '{other}' for {context} (expected human, mouse, yeast, or ecoli)"
        )),
    }
}

fn parse_translation_speed_mark_for_shell(
    raw: &str,
    context: &str,
) -> Result<TranslationSpeedMark, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "fast" => Ok(TranslationSpeedMark::Fast),
        "slow" => Ok(TranslationSpeedMark::Slow),
        other => Err(format!(
            "Unknown --speed-mark '{other}' for {context} (expected fast or slow)"
        )),
    }
}

fn parse_protein_to_dna_handoff_ranking_goal(
    raw: &str,
    context: &str,
) -> Result<ProteinToDnaHandoffRankingGoal, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "balanced" | "balanced_provenance" | "balanced-provenance" => {
            Ok(ProteinToDnaHandoffRankingGoal::BalancedProvenance)
        }
        "native" | "native_fidelity" | "native-fidelity" => {
            Ok(ProteinToDnaHandoffRankingGoal::NativeFidelity)
        }
        "expression" | "expression_optimized" | "expression-optimized" | "optimized" => {
            Ok(ProteinToDnaHandoffRankingGoal::ExpressionOptimized)
        }
        other => Err(format!(
            "Unsupported ranking goal '{other}' for {context}; expected balanced_provenance, native_fidelity, or expression_optimized"
        )),
    }
}

pub(super) fn parse_reverse_translate_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "reverse-translate requires a subcommand: run, list-reports, show-report, export-report"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "run" => {
            if tokens.len() < 3 {
                return Err(
                    "reverse-translate run requires PROTEIN_SEQ_ID [--output-id ID] [--speed-profile human|mouse|yeast|ecoli] [--speed-mark fast|slow] [--translation-table N] [--target-anneal-tm-c F] [--anneal-window-bp N]"
                        .to_string(),
                );
            }
            let seq_id = tokens[2].trim().to_string();
            if seq_id.is_empty() {
                return Err("reverse-translate run PROTEIN_SEQ_ID must not be empty".to_string());
            }
            let mut output_id: Option<String> = None;
            let mut speed_profile: Option<TranslationSpeedProfile> = None;
            let mut speed_mark: Option<TranslationSpeedMark> = None;
            let mut translation_table: Option<usize> = None;
            let mut target_anneal_tm_c: Option<f64> = None;
            let mut anneal_window_bp: Option<usize> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--output-id" | "--output" => {
                        let flag = tokens[idx].clone();
                        let raw =
                            parse_option_path(tokens, &mut idx, &flag, "reverse-translate run")?;
                        let trimmed = raw.trim();
                        output_id = if trimmed.is_empty() {
                            None
                        } else {
                            Some(trimmed.to_string())
                        };
                    }
                    "--speed-profile" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--speed-profile",
                            "reverse-translate run",
                        )?;
                        speed_profile = Some(parse_translation_speed_profile_for_shell(
                            &raw,
                            "reverse-translate run",
                        )?);
                    }
                    "--speed-mark" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--speed-mark",
                            "reverse-translate run",
                        )?;
                        speed_mark = Some(parse_translation_speed_mark_for_shell(
                            &raw,
                            "reverse-translate run",
                        )?);
                    }
                    "--translation-table" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--translation-table",
                            "reverse-translate run",
                        )?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --translation-table value '{raw}' for reverse-translate run: {e}"
                            )
                        })?;
                        if parsed == 0 {
                            return Err(
                                "--translation-table for reverse-translate run must be >= 1"
                                    .to_string(),
                            );
                        }
                        translation_table = Some(parsed);
                    }
                    "--target-anneal-tm-c" | "--target-anneal-tm" => {
                        let flag = tokens[idx].clone();
                        let raw =
                            parse_option_path(tokens, &mut idx, &flag, "reverse-translate run")?;
                        target_anneal_tm_c = Some(raw.parse::<f64>().map_err(|e| {
                            format!("Invalid {flag} value '{raw}' for reverse-translate run: {e}")
                        })?);
                    }
                    "--anneal-window-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--anneal-window-bp",
                            "reverse-translate run",
                        )?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --anneal-window-bp value '{raw}' for reverse-translate run: {e}"
                            )
                        })?;
                        if parsed == 0 {
                            return Err(
                                "--anneal-window-bp for reverse-translate run must be >= 1"
                                    .to_string(),
                            );
                        }
                        anneal_window_bp = Some(parsed);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for reverse-translate run"
                        ));
                    }
                }
            }
            Ok(ShellCommand::ReverseTranslateRun {
                seq_id,
                output_id,
                speed_profile,
                speed_mark,
                translation_table,
                target_anneal_tm_c,
                anneal_window_bp,
            })
        }
        "list-reports" => {
            let mut protein_seq_id: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--protein" | "--protein-seq-id" | "--seq-id" => {
                        let flag = tokens[idx].clone();
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            &flag,
                            "reverse-translate list-reports",
                        )?;
                        let trimmed = raw.trim();
                        protein_seq_id = if trimmed.is_empty() {
                            None
                        } else {
                            Some(trimmed.to_string())
                        };
                    }
                    other => {
                        if protein_seq_id.is_none() && !other.starts_with('-') {
                            protein_seq_id = Some(other.to_string());
                            idx += 1;
                        } else {
                            return Err(format!(
                                "Unknown option '{other}' for reverse-translate list-reports"
                            ));
                        }
                    }
                }
            }
            Ok(ShellCommand::ReverseTranslateListReports { protein_seq_id })
        }
        "show-report" => {
            if tokens.len() != 3 {
                return Err("reverse-translate show-report requires REPORT_ID".to_string());
            }
            Ok(ShellCommand::ReverseTranslateShowReport {
                report_id: tokens[2].trim().to_string(),
            })
        }
        "export-report" => {
            if tokens.len() != 4 {
                return Err(
                    "reverse-translate export-report requires REPORT_ID OUTPUT.json".to_string(),
                );
            }
            Ok(ShellCommand::ReverseTranslateExportReport {
                report_id: tokens[2].trim().to_string(),
                path: tokens[3].clone(),
            })
        }
        other => Err(format!(
            "Unknown reverse-translate subcommand '{other}' (expected run, list-reports, show-report, export-report)"
        )),
    }
}

pub(super) fn parse_construct_reasoning_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "construct-reasoning requires a subcommand: build-protein-dna-handoff, list-graphs, show-graph, set-annotation-status, export-graph"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "build-protein-dna-handoff" | "build-protein-handoff" => {
            if tokens.len() < 4 {
                return Err(
                    "construct-reasoning build-protein-dna-handoff requires SEQ_ID PROTEIN_SEQ_ID [--transcript TRANSCRIPT_ID] [--projection-id ID] [--ensembl-entry ID] [--feature-query TEXT] [--ranking-goal GOAL] [--speed-profile PROFILE] [--speed-mark MARK] [--translation-table N] [--target-anneal-tm-c N] [--anneal-window-bp N] [--objective-id ID] [--graph-id ID]"
                        .to_string(),
                );
            }
            let seq_id = tokens[2].trim().to_string();
            let protein_seq_id = tokens[3].trim().to_string();
            if seq_id.is_empty() || protein_seq_id.is_empty() {
                return Err(
                    "construct-reasoning build-protein-dna-handoff requires non-empty SEQ_ID and PROTEIN_SEQ_ID"
                        .to_string(),
                );
            }
            let mut transcript_filter: Option<String> = None;
            let mut projection_id: Option<String> = None;
            let mut ensembl_entry_id: Option<String> = None;
            let mut feature_query: Option<String> = None;
            let mut ranking_goal = ProteinToDnaHandoffRankingGoal::BalancedProvenance;
            let mut speed_profile: Option<TranslationSpeedProfile> = None;
            let mut speed_mark: Option<TranslationSpeedMark> = None;
            let mut translation_table: Option<usize> = None;
            let mut target_anneal_tm_c: Option<f64> = None;
            let mut anneal_window_bp: Option<usize> = None;
            let mut objective_id: Option<String> = None;
            let mut graph_id: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--transcript" | "--transcript-filter" => {
                        let flag = tokens[idx].clone();
                        transcript_filter = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            &flag,
                            "construct-reasoning build-protein-dna-handoff",
                        )?);
                    }
                    "--projection-id" | "--projection" => {
                        let flag = tokens[idx].clone();
                        projection_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            &flag,
                            "construct-reasoning build-protein-dna-handoff",
                        )?);
                    }
                    "--ensembl-entry" | "--ensembl-entry-id" => {
                        let flag = tokens[idx].clone();
                        ensembl_entry_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            &flag,
                            "construct-reasoning build-protein-dna-handoff",
                        )?);
                    }
                    "--feature-query" | "--feature" => {
                        let flag = tokens[idx].clone();
                        feature_query = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            &flag,
                            "construct-reasoning build-protein-dna-handoff",
                        )?);
                    }
                    "--ranking-goal" | "--goal" => {
                        let flag = tokens[idx].clone();
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            &flag,
                            "construct-reasoning build-protein-dna-handoff",
                        )?;
                        ranking_goal = parse_protein_to_dna_handoff_ranking_goal(
                            &raw,
                            "construct-reasoning build-protein-dna-handoff",
                        )?;
                    }
                    "--speed-profile" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--speed-profile",
                            "construct-reasoning build-protein-dna-handoff",
                        )?;
                        speed_profile = Some(parse_translation_speed_profile_for_shell(
                            &raw,
                            "construct-reasoning build-protein-dna-handoff",
                        )?);
                    }
                    "--speed-mark" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--speed-mark",
                            "construct-reasoning build-protein-dna-handoff",
                        )?;
                        speed_mark = Some(parse_translation_speed_mark_for_shell(
                            &raw,
                            "construct-reasoning build-protein-dna-handoff",
                        )?);
                    }
                    "--translation-table" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--translation-table",
                            "construct-reasoning build-protein-dna-handoff",
                        )?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --translation-table value '{raw}' for construct-reasoning build-protein-dna-handoff: {e}"
                            )
                        })?;
                        if parsed == 0 {
                            return Err(
                                "--translation-table for construct-reasoning build-protein-dna-handoff must be >= 1"
                                    .to_string(),
                            );
                        }
                        translation_table = Some(parsed);
                    }
                    "--target-anneal-tm-c" | "--target-anneal-tm" => {
                        let flag = tokens[idx].clone();
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            &flag,
                            "construct-reasoning build-protein-dna-handoff",
                        )?;
                        target_anneal_tm_c = Some(raw.parse::<f64>().map_err(|e| {
                            format!(
                                "Invalid {flag} value '{raw}' for construct-reasoning build-protein-dna-handoff: {e}"
                            )
                        })?);
                    }
                    "--anneal-window-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--anneal-window-bp",
                            "construct-reasoning build-protein-dna-handoff",
                        )?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --anneal-window-bp value '{raw}' for construct-reasoning build-protein-dna-handoff: {e}"
                            )
                        })?;
                        if parsed == 0 {
                            return Err(
                                "--anneal-window-bp for construct-reasoning build-protein-dna-handoff must be >= 1"
                                    .to_string(),
                            );
                        }
                        anneal_window_bp = Some(parsed);
                    }
                    "--objective-id" => {
                        objective_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--objective-id",
                            "construct-reasoning build-protein-dna-handoff",
                        )?);
                    }
                    "--graph-id" => {
                        graph_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--graph-id",
                            "construct-reasoning build-protein-dna-handoff",
                        )?);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for construct-reasoning build-protein-dna-handoff"
                        ));
                    }
                }
            }
            Ok(ShellCommand::ConstructReasoningBuildProteinDnaHandoff {
                seq_id,
                protein_seq_id,
                transcript_filter,
                projection_id,
                ensembl_entry_id,
                feature_query,
                ranking_goal,
                speed_profile,
                speed_mark,
                translation_table,
                target_anneal_tm_c,
                anneal_window_bp,
                objective_id,
                graph_id,
            })
        }
        "list-graphs" => {
            let mut seq_id: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--seq-id" => {
                        seq_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--seq-id",
                            "construct-reasoning list-graphs",
                        )?);
                    }
                    other => {
                        if seq_id.is_none() && !other.starts_with('-') {
                            seq_id = Some(other.to_string());
                            idx += 1;
                        } else {
                            return Err(format!(
                                "Unknown option '{other}' for construct-reasoning list-graphs"
                            ));
                        }
                    }
                }
            }
            Ok(ShellCommand::ConstructReasoningListGraphs { seq_id })
        }
        "show-graph" => {
            if tokens.len() != 3 {
                return Err("construct-reasoning show-graph requires GRAPH_ID".to_string());
            }
            Ok(ShellCommand::ConstructReasoningShowGraph {
                graph_id: tokens[2].trim().to_string(),
            })
        }
        "set-annotation-status" => {
            if tokens.len() != 5 {
                return Err(
                    "construct-reasoning set-annotation-status requires GRAPH_ID ANNOTATION_ID STATUS"
                        .to_string(),
                );
            }
            let graph_id = tokens[2].trim().to_string();
            let annotation_id = tokens[3].trim().to_string();
            if graph_id.is_empty() || annotation_id.is_empty() {
                return Err(
                    "construct-reasoning set-annotation-status requires non-empty GRAPH_ID and ANNOTATION_ID"
                        .to_string(),
                );
            }
            let editable_status = match tokens[4].trim().to_ascii_lowercase().as_str() {
                "draft" => EditableStatus::Draft,
                "accepted" => EditableStatus::Accepted,
                "rejected" => EditableStatus::Rejected,
                "locked" => EditableStatus::Locked,
                other => {
                    return Err(format!(
                        "Unsupported construct-reasoning annotation status '{other}' (expected draft|accepted|rejected|locked)"
                    ));
                }
            };
            Ok(ShellCommand::ConstructReasoningSetAnnotationStatus {
                graph_id,
                annotation_id,
                editable_status,
            })
        }
        "write-annotation" => {
            if tokens.len() != 4 {
                return Err(
                    "construct-reasoning write-annotation requires GRAPH_ID ANNOTATION_ID"
                        .to_string(),
                );
            }
            let graph_id = tokens[2].trim().to_string();
            let annotation_id = tokens[3].trim().to_string();
            if graph_id.is_empty() || annotation_id.is_empty() {
                return Err(
                    "construct-reasoning write-annotation requires non-empty GRAPH_ID and ANNOTATION_ID"
                        .to_string(),
                );
            }
            Ok(ShellCommand::ConstructReasoningWriteAnnotation {
                graph_id,
                annotation_id,
            })
        }
        "export-graph" => {
            if tokens.len() != 4 {
                return Err(
                    "construct-reasoning export-graph requires GRAPH_ID OUTPUT.json".to_string(),
                );
            }
            Ok(ShellCommand::ConstructReasoningExportGraph {
                graph_id: tokens[2].trim().to_string(),
                path: tokens[3].clone(),
            })
        }
        other => Err(format!(
            "Unknown construct-reasoning subcommand '{other}' (expected build-protein-dna-handoff, list-graphs, show-graph, set-annotation-status, write-annotation, export-graph)"
        )),
    }
}

pub(super) fn parse_seq_primer_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("seq-primer requires a subcommand: suggest".to_string());
    }
    match tokens[1].as_str() {
        "suggest" => {
            if tokens.len() < 4 {
                return Err(
                    "seq-primer suggest requires EXPECTED_SEQ_ID plus optional --primers/--primer values and/or --confirmation-report REPORT_ID [--min-3prime-anneal-bp N] [--predicted-read-length-bp N]"
                        .to_string(),
                );
            }
            let expected_seq_id = tokens[2].trim().to_string();
            if expected_seq_id.is_empty() {
                return Err("seq-primer suggest EXPECTED_SEQ_ID must not be empty".to_string());
            }
            let mut primer_seq_ids: Vec<String> = vec![];
            let mut confirmation_report_id: Option<String> = None;
            let mut min_3prime_anneal_bp = 18usize;
            let mut predicted_read_length_bp = 800usize;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--primers" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--primers", "seq-primer suggest")?;
                        let parsed = raw
                            .split(',')
                            .map(str::trim)
                            .filter(|value| !value.is_empty())
                            .map(|value| value.to_string())
                            .collect::<Vec<_>>();
                        if parsed.is_empty() {
                            return Err(
                                "--primers for seq-primer suggest must include at least one ID"
                                    .to_string(),
                            );
                        }
                        primer_seq_ids.extend(parsed);
                    }
                    "--primer" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--primer", "seq-primer suggest")?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            return Err(
                                "--primer for seq-primer suggest must not be empty".to_string()
                            );
                        }
                        primer_seq_ids.push(trimmed.to_string());
                    }
                    "--confirmation-report" | "--report-id" => {
                        let flag = tokens[idx].clone();
                        let raw = parse_option_path(tokens, &mut idx, &flag, "seq-primer suggest")?;
                        let trimmed = raw.trim();
                        confirmation_report_id = if trimmed.is_empty() {
                            None
                        } else {
                            Some(trimmed.to_string())
                        };
                    }
                    "--min-3prime-anneal-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-3prime-anneal-bp",
                            "seq-primer suggest",
                        )?;
                        min_3prime_anneal_bp = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --min-3prime-anneal-bp value '{raw}' for seq-primer suggest: {e}"
                            )
                        })?;
                    }
                    "--predicted-read-length-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--predicted-read-length-bp",
                            "seq-primer suggest",
                        )?;
                        predicted_read_length_bp = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --predicted-read-length-bp value '{raw}' for seq-primer suggest: {e}"
                            )
                        })?;
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for seq-primer suggest"));
                    }
                }
            }
            if primer_seq_ids.is_empty() && confirmation_report_id.is_none() {
                return Err(
                    "seq-primer suggest requires at least one primer via --primers/--primer or a saved report via --confirmation-report"
                        .to_string(),
                );
            }
            Ok(ShellCommand::SeqPrimerSuggest {
                expected_seq_id,
                primer_seq_ids,
                confirmation_report_id,
                min_3prime_anneal_bp,
                predicted_read_length_bp,
            })
        }
        other => Err(format!(
            "Unknown seq-primer subcommand '{other}' (expected suggest)"
        )),
    }
}

pub(super) fn parse_rna_reads_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "rna-reads requires a subcommand: interpret, align-report, list-reports, show-report, summarize-gene-support, inspect-gene-support, inspect-alignments, inspect-concatemers, export-report, export-hits-fasta, export-sample-sheet, export-paths-tsv, export-abundance-tsv, export-score-density-svg, export-alignments-tsv, export-alignment-dotplot-svg"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "interpret" => {
            if tokens.len() < 5 {
                return Err(
                    "rna-reads interpret requires SEQ_ID FEATURE_ID INPUT.fa[.gz] [--report-id ID] [--report-mode full|seed_passed_only] [--checkpoint-path PATH] [--checkpoint-every-reads N] [--resume-from-checkpoint|--no-resume-from-checkpoint] [--profile PROFILE] [--format fasta] [--scope SCOPE] [--origin-mode single_gene|multi_gene_sparse] [--target-gene GENE_ID]... [--roi-seed-capture|--no-roi-seed-capture] [--kmer-len N] [--seed-stride-bp N] [--min-seed-hit-fraction F] [--min-weighted-seed-hit-fraction F] [--min-unique-matched-kmers N] [--min-chain-consistency-fraction F] [--max-median-transcript-gap F] [--min-confirmed-transitions N] [--min-transition-support-fraction F] [--cdna-poly-t-flip|--no-cdna-poly-t-flip] [--poly-t-prefix-min-bp N] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]"
                        .to_string(),
                );
            }
            let seq_id = tokens[2].trim().to_string();
            if seq_id.is_empty() {
                return Err("rna-reads interpret SEQ_ID must not be empty".to_string());
            }
            let seed_feature_id = tokens[3].parse::<usize>().map_err(|e| {
                format!(
                    "Invalid FEATURE_ID '{}' for rna-reads interpret: {e}",
                    tokens[3]
                )
            })?;
            let input_path = tokens[4].trim().to_string();
            if input_path.is_empty() {
                return Err("rna-reads interpret INPUT.fa[.gz] must not be empty".to_string());
            }
            let mut profile = RnaReadInterpretationProfile::NanoporeCdnaV1;
            let mut input_format = RnaReadInputFormat::Fasta;
            let mut scope = SplicingScopePreset::AllOverlappingBothStrands;
            let mut origin_mode = RnaReadOriginMode::SingleGene;
            let mut target_gene_ids: Vec<String> = vec![];
            let mut roi_seed_capture_enabled = false;
            let mut seed_filter = RnaReadSeedFilterConfig::default();
            let mut align_config = RnaReadAlignConfig::default();
            let mut report_id: Option<String> = None;
            let mut report_mode = RnaReadReportMode::Full;
            let mut checkpoint_path: Option<String> = None;
            let mut checkpoint_every_reads = 10_000usize;
            let mut resume_from_checkpoint = false;
            let mut idx = 5usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--report-id" => {
                        report_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--report-id",
                            "rna-reads interpret",
                        )?);
                    }
                    "--profile" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--profile",
                            "rna-reads interpret",
                        )?;
                        profile = parse_rna_read_profile(&raw)?;
                    }
                    "--report-mode" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--report-mode",
                            "rna-reads interpret",
                        )?;
                        report_mode = parse_rna_read_report_mode(&raw)?;
                    }
                    "--checkpoint-path" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--checkpoint-path",
                            "rna-reads interpret",
                        )?;
                        checkpoint_path = Some(raw);
                    }
                    "--checkpoint-every-reads" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--checkpoint-every-reads",
                            "rna-reads interpret",
                        )?;
                        checkpoint_every_reads = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --checkpoint-every-reads value '{raw}' for rna-reads interpret: {e}"
                            )
                        })?;
                    }
                    "--resume-from-checkpoint" => {
                        resume_from_checkpoint = true;
                        idx += 1;
                    }
                    "--no-resume-from-checkpoint" => {
                        resume_from_checkpoint = false;
                        idx += 1;
                    }
                    "--format" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--format", "rna-reads interpret")?;
                        input_format = parse_rna_read_input_format(&raw)?;
                    }
                    "--scope" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--scope", "rna-reads interpret")?;
                        scope = parse_splicing_scope_preset(&raw)?;
                    }
                    "--origin-mode" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--origin-mode",
                            "rna-reads interpret",
                        )?;
                        origin_mode = parse_rna_read_origin_mode(&raw)?;
                    }
                    "--target-gene" | "--target-gene-id" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--target-gene",
                            "rna-reads interpret",
                        )?;
                        let gene_id = raw.trim();
                        if gene_id.is_empty() {
                            return Err(
                                "--target-gene requires a non-empty gene identifier".to_string()
                            );
                        }
                        target_gene_ids.push(gene_id.to_string());
                    }
                    "--roi-seed-capture" => {
                        roi_seed_capture_enabled = true;
                        idx += 1;
                    }
                    "--no-roi-seed-capture" => {
                        roi_seed_capture_enabled = false;
                        idx += 1;
                    }
                    "--kmer-len" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--kmer-len",
                            "rna-reads interpret",
                        )?;
                        seed_filter.kmer_len = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --kmer-len value '{raw}' for rna-reads interpret: {e}")
                        })?;
                    }
                    "--seed-stride-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--seed-stride-bp",
                            "rna-reads interpret",
                        )?;
                        seed_filter.seed_stride_bp = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --seed-stride-bp value '{raw}' for rna-reads interpret: {e}"
                            )
                        })?;
                    }
                    "--min-seed-hit-fraction" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-seed-hit-fraction",
                            "rna-reads interpret",
                        )?;
                        seed_filter.min_seed_hit_fraction =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --min-seed-hit-fraction value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--min-weighted-seed-hit-fraction" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-weighted-seed-hit-fraction",
                            "rna-reads interpret",
                        )?;
                        seed_filter.min_weighted_seed_hit_fraction =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --min-weighted-seed-hit-fraction value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--min-unique-matched-kmers" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-unique-matched-kmers",
                            "rna-reads interpret",
                        )?;
                        seed_filter.min_unique_matched_kmers =
                            raw.parse::<usize>().map_err(|e| {
                                format!(
                                    "Invalid --min-unique-matched-kmers value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--max-median-transcript-gap" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-median-transcript-gap",
                            "rna-reads interpret",
                        )?;
                        seed_filter.max_median_transcript_gap =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --max-median-transcript-gap value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--min-chain-consistency-fraction" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-chain-consistency-fraction",
                            "rna-reads interpret",
                        )?;
                        seed_filter.min_chain_consistency_fraction =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --min-chain-consistency-fraction value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--min-confirmed-transitions" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-confirmed-transitions",
                            "rna-reads interpret",
                        )?;
                        seed_filter.min_confirmed_exon_transitions =
                            raw.parse::<usize>().map_err(|e| {
                                format!(
                                    "Invalid --min-confirmed-transitions value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--min-transition-support-fraction" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-transition-support-fraction",
                            "rna-reads interpret",
                        )?;
                        seed_filter.min_transition_support_fraction =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --min-transition-support-fraction value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--cdna-poly-t-flip" => {
                        seed_filter.cdna_poly_t_flip_enabled = true;
                        idx += 1;
                    }
                    "--no-cdna-poly-t-flip" => {
                        seed_filter.cdna_poly_t_flip_enabled = false;
                        idx += 1;
                    }
                    "--poly-t-prefix-min-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--poly-t-prefix-min-bp",
                            "rna-reads interpret",
                        )?;
                        seed_filter.poly_t_prefix_min_bp = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --poly-t-prefix-min-bp value '{raw}' for rna-reads interpret: {e}"
                            )
                        })?;
                    }
                    "--align-band-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--align-band-bp",
                            "rna-reads interpret",
                        )?;
                        align_config.band_width_bp = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --align-band-bp value '{raw}' for rna-reads interpret: {e}"
                            )
                        })?;
                    }
                    "--align-min-identity" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--align-min-identity",
                            "rna-reads interpret",
                        )?;
                        align_config.min_identity_fraction = raw.parse::<f64>().map_err(|e| {
                            format!(
                                "Invalid --align-min-identity value '{raw}' for rna-reads interpret: {e}"
                            )
                        })?;
                    }
                    "--max-secondary-mappings" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-secondary-mappings",
                            "rna-reads interpret",
                        )?;
                        align_config.max_secondary_mappings =
                            raw.parse::<usize>().map_err(|e| {
                                format!(
                                    "Invalid --max-secondary-mappings value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for rna-reads interpret"));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsInterpret {
                seq_id,
                seed_feature_id,
                input_path,
                profile,
                input_format,
                scope,
                origin_mode,
                target_gene_ids,
                roi_seed_capture_enabled,
                seed_filter,
                align_config,
                report_id,
                report_mode,
                checkpoint_path,
                checkpoint_every_reads,
                resume_from_checkpoint,
            })
        }
        "align-report" => {
            if tokens.len() < 3 {
                return Err(
                    "rna-reads align-report requires REPORT_ID [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].trim().to_string();
            if report_id.is_empty() {
                return Err("rna-reads align-report REPORT_ID must not be empty".to_string());
            }
            let mut selection = RnaReadHitSelection::SeedPassed;
            let mut align_config_override: Option<RnaReadAlignConfig> = None;
            let mut selected_record_indices: Vec<usize> = vec![];
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--selection" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--selection",
                            "rna-reads align-report",
                        )?;
                        selection = parse_rna_read_hit_selection(&raw)?;
                    }
                    "--record-indices" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--record-indices",
                            "rna-reads align-report",
                        )?;
                        selected_record_indices = parse_rna_read_record_indices(&raw)?;
                    }
                    "--align-band-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--align-band-bp",
                            "rna-reads align-report",
                        )?;
                        let cfg =
                            align_config_override.get_or_insert_with(RnaReadAlignConfig::default);
                        cfg.band_width_bp = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --align-band-bp value '{raw}' for rna-reads align-report: {e}"
                            )
                        })?;
                    }
                    "--align-min-identity" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--align-min-identity",
                            "rna-reads align-report",
                        )?;
                        let cfg =
                            align_config_override.get_or_insert_with(RnaReadAlignConfig::default);
                        cfg.min_identity_fraction = raw.parse::<f64>().map_err(|e| {
                            format!(
                                "Invalid --align-min-identity value '{raw}' for rna-reads align-report: {e}"
                            )
                        })?;
                    }
                    "--max-secondary-mappings" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-secondary-mappings",
                            "rna-reads align-report",
                        )?;
                        let cfg =
                            align_config_override.get_or_insert_with(RnaReadAlignConfig::default);
                        cfg.max_secondary_mappings = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --max-secondary-mappings value '{raw}' for rna-reads align-report: {e}"
                            )
                        })?;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads align-report"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsAlignReport {
                report_id,
                selection,
                align_config_override,
                selected_record_indices,
            })
        }
        "list-reports" => {
            if tokens.len() > 3 {
                return Err(
                    "rna-reads list-reports expects at most one optional SEQ_ID".to_string()
                );
            }
            let seq_id = if tokens.len() == 3 {
                let value = tokens[2].trim();
                if value.is_empty() {
                    None
                } else {
                    Some(value.to_string())
                }
            } else {
                None
            };
            Ok(ShellCommand::RnaReadsListReports { seq_id })
        }
        "show-report" => {
            if tokens.len() != 3 {
                return Err("rna-reads show-report requires REPORT_ID".to_string());
            }
            Ok(ShellCommand::RnaReadsShowReport {
                report_id: tokens[2].clone(),
            })
        }
        "summarize-gene-support" => {
            if tokens.len() < 3 {
                return Err(
                    "rna-reads summarize-gene-support requires REPORT_ID --gene GENE [--gene GENE ...] [--record-indices i,j,k] [--complete-rule near|strict|exact] [--output PATH]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].trim().to_string();
            if report_id.is_empty() {
                return Err(
                    "rna-reads summarize-gene-support REPORT_ID must not be empty".to_string(),
                );
            }
            let mut gene_ids = Vec::<String>::new();
            let mut selected_record_indices = Vec::<usize>::new();
            let mut complete_rule = RnaReadGeneSupportCompleteRule::Near;
            let mut output_path: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--gene" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--gene",
                            "rna-reads summarize-gene-support",
                        )?;
                        gene_ids.push(raw);
                    }
                    "--record-indices" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--record-indices",
                            "rna-reads summarize-gene-support",
                        )?;
                        selected_record_indices = parse_rna_read_record_indices(&raw)?;
                    }
                    "--complete-rule" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--complete-rule",
                            "rna-reads summarize-gene-support",
                        )?;
                        complete_rule = parse_rna_read_gene_support_complete_rule(&raw)?;
                    }
                    "--output" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--output",
                            "rna-reads summarize-gene-support",
                        )?;
                        output_path = Some(raw);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads summarize-gene-support"
                        ));
                    }
                }
            }
            if gene_ids.is_empty() {
                return Err(
                    "rna-reads summarize-gene-support requires at least one --gene GENE"
                        .to_string(),
                );
            }
            Ok(ShellCommand::RnaReadsSummarizeGeneSupport {
                report_id,
                gene_ids,
                selected_record_indices,
                complete_rule,
                output_path,
            })
        }
        "inspect-gene-support" => {
            if tokens.len() < 3 {
                return Err(
                    "rna-reads inspect-gene-support requires REPORT_ID --gene GENE [--gene GENE ...] [--record-indices i,j,k] [--complete-rule near|strict|exact] [--cohort all|accepted|fragment|complete|rejected] [--output PATH]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].trim().to_string();
            if report_id.is_empty() {
                return Err(
                    "rna-reads inspect-gene-support REPORT_ID must not be empty".to_string()
                );
            }
            let mut gene_ids = Vec::<String>::new();
            let mut selected_record_indices = Vec::<usize>::new();
            let mut complete_rule = RnaReadGeneSupportCompleteRule::Near;
            let mut cohort_filter = RnaReadGeneSupportAuditCohortFilter::All;
            let mut output_path: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--gene" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--gene",
                            "rna-reads inspect-gene-support",
                        )?;
                        gene_ids.push(raw);
                    }
                    "--record-indices" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--record-indices",
                            "rna-reads inspect-gene-support",
                        )?;
                        selected_record_indices = parse_rna_read_record_indices(&raw)?;
                    }
                    "--complete-rule" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--complete-rule",
                            "rna-reads inspect-gene-support",
                        )?;
                        complete_rule = parse_rna_read_gene_support_complete_rule(&raw)?;
                    }
                    "--cohort" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--cohort",
                            "rna-reads inspect-gene-support",
                        )?;
                        cohort_filter = parse_rna_read_gene_support_audit_cohort_filter(&raw)?;
                    }
                    "--output" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--output",
                            "rna-reads inspect-gene-support",
                        )?;
                        output_path = Some(raw);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads inspect-gene-support"
                        ));
                    }
                }
            }
            if gene_ids.is_empty() {
                return Err(
                    "rna-reads inspect-gene-support requires at least one --gene GENE".to_string(),
                );
            }
            Ok(ShellCommand::RnaReadsInspectGeneSupport {
                report_id,
                gene_ids,
                selected_record_indices,
                complete_rule,
                cohort_filter,
                output_path,
            })
        }
        "inspect-alignments" => {
            if tokens.len() < 3 {
                return Err(
                    "rna-reads inspect-alignments requires REPORT_ID [--selection all|seed_passed|aligned] [--limit N] [--effect-filter all_aligned|confirmed_only|disagreement_only|reassigned_only|no_phase1_only|selected_only] [--sort rank|identity|coverage|score] [--search TEXT] [--record-indices i,j,k] [--score-bin-index N] [--score-bin-count M]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].trim().to_string();
            if report_id.is_empty() {
                return Err("rna-reads inspect-alignments REPORT_ID must not be empty".to_string());
            }
            let mut selection = RnaReadHitSelection::Aligned;
            let mut limit = 50usize;
            let mut effect_filter = RnaReadAlignmentInspectionEffectFilter::AllAligned;
            let mut sort_key = RnaReadAlignmentInspectionSortKey::Rank;
            let mut search = String::new();
            let mut selected_record_indices: Vec<usize> = vec![];
            let mut score_density_variant = RnaReadScoreDensityVariant::AllScored;
            let mut score_bin_index: Option<usize> = None;
            let mut score_bin_count = 0usize;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--selection" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--selection",
                            "rna-reads inspect-alignments",
                        )?;
                        selection = parse_rna_read_hit_selection(&raw)?;
                    }
                    "--limit" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--limit",
                            "rna-reads inspect-alignments",
                        )?;
                        limit = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --limit value '{raw}' for rna-reads inspect-alignments: {e}"
                            )
                        })?;
                    }
                    "--effect-filter" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--effect-filter",
                            "rna-reads inspect-alignments",
                        )?;
                        effect_filter = parse_rna_read_alignment_effect_filter(&raw)?;
                    }
                    "--sort" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--sort",
                            "rna-reads inspect-alignments",
                        )?;
                        sort_key = parse_rna_read_alignment_sort_key(&raw)?;
                    }
                    "--search" => {
                        search = parse_option_path(
                            tokens,
                            &mut idx,
                            "--search",
                            "rna-reads inspect-alignments",
                        )?;
                    }
                    "--score-bin-variant" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--score-bin-variant",
                            "rna-reads inspect-alignments",
                        )?;
                        score_density_variant = parse_rna_read_score_density_variant(&raw)?;
                    }
                    "--record-indices" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--record-indices",
                            "rna-reads inspect-alignments",
                        )?;
                        selected_record_indices = parse_rna_read_record_indices(&raw)?;
                    }
                    "--score-bin-index" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--score-bin-index",
                            "rna-reads inspect-alignments",
                        )?;
                        score_bin_index = Some(raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --score-bin-index value '{raw}' for rna-reads inspect-alignments: {e}"
                            )
                        })?);
                    }
                    "--score-bin-count" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--score-bin-count",
                            "rna-reads inspect-alignments",
                        )?;
                        score_bin_count = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --score-bin-count value '{raw}' for rna-reads inspect-alignments: {e}"
                            )
                        })?;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads inspect-alignments"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsInspectAlignments {
                report_id,
                selection,
                limit,
                effect_filter,
                sort_key,
                search,
                selected_record_indices,
                score_density_variant,
                score_bin_index,
                score_bin_count,
            })
        }
        "inspect-concatemers" => {
            if tokens.len() < 3 {
                return Err(
                    "rna-reads inspect-concatemers requires REPORT_ID [--selection all|seed_passed|aligned] [--limit N] [--record-indices i,j,k] [--internal-homopolymer-min-bp N] [--end-margin-bp N] [--max-primary-query-cov F] [--min-secondary-identity F] [--max-secondary-query-overlap F] [--adapter-fasta PATH] [--adapter-min-match-bp N] [--fragment-min-bp N] [--fragment-max-parts N] [--fragment-min-identity F] [--fragment-min-query-cov F] [--transcript-fasta PATH]... [--transcript-index PATH]..."
                        .to_string(),
                );
            }
            let report_id = tokens[2].trim().to_string();
            if report_id.is_empty() {
                return Err("rna-reads inspect-concatemers REPORT_ID must not be empty".to_string());
            }
            let mut selection = RnaReadHitSelection::All;
            let mut limit = 50usize;
            let mut selected_record_indices: Vec<usize> = vec![];
            let mut settings = RnaReadConcatemerInspectionSettings::default();
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--selection" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--selection",
                            "rna-reads inspect-concatemers",
                        )?;
                        selection = parse_rna_read_hit_selection(&raw)?;
                    }
                    "--limit" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--limit",
                            "rna-reads inspect-concatemers",
                        )?;
                        limit = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --limit value '{raw}' for rna-reads inspect-concatemers: {e}"
                            )
                        })?;
                    }
                    "--record-indices" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--record-indices",
                            "rna-reads inspect-concatemers",
                        )?;
                        selected_record_indices = parse_rna_read_record_indices(&raw)?;
                    }
                    "--internal-homopolymer-min-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--internal-homopolymer-min-bp",
                            "rna-reads inspect-concatemers",
                        )?;
                        settings.internal_homopolymer_min_bp =
                            raw.parse::<usize>().map_err(|e| {
                                format!(
                                    "Invalid --internal-homopolymer-min-bp value '{raw}' for rna-reads inspect-concatemers: {e}"
                                )
                            })?;
                    }
                    "--end-margin-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--end-margin-bp",
                            "rna-reads inspect-concatemers",
                        )?;
                        settings.end_margin_bp = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --end-margin-bp value '{raw}' for rna-reads inspect-concatemers: {e}"
                            )
                        })?;
                    }
                    "--max-primary-query-cov" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-primary-query-cov",
                            "rna-reads inspect-concatemers",
                        )?;
                        settings.max_primary_query_coverage_fraction =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --max-primary-query-cov value '{raw}' for rna-reads inspect-concatemers: {e}"
                                )
                            })?;
                    }
                    "--min-secondary-identity" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-secondary-identity",
                            "rna-reads inspect-concatemers",
                        )?;
                        settings.min_secondary_identity_fraction =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --min-secondary-identity value '{raw}' for rna-reads inspect-concatemers: {e}"
                                )
                            })?;
                    }
                    "--max-secondary-query-overlap" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-secondary-query-overlap",
                            "rna-reads inspect-concatemers",
                        )?;
                        settings.max_secondary_query_overlap_fraction =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --max-secondary-query-overlap value '{raw}' for rna-reads inspect-concatemers: {e}"
                                )
                            })?;
                    }
                    "--adapter-fasta" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--adapter-fasta",
                            "rna-reads inspect-concatemers",
                        )?;
                        settings.adapter_fasta_path = Some(raw);
                    }
                    "--adapter-min-match-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--adapter-min-match-bp",
                            "rna-reads inspect-concatemers",
                        )?;
                        settings.adapter_min_match_bp = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --adapter-min-match-bp value '{raw}' for rna-reads inspect-concatemers: {e}"
                            )
                        })?;
                    }
                    "--fragment-min-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--fragment-min-bp",
                            "rna-reads inspect-concatemers",
                        )?;
                        settings.fragment_min_bp = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --fragment-min-bp value '{raw}' for rna-reads inspect-concatemers: {e}"
                            )
                        })?;
                    }
                    "--fragment-max-parts" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--fragment-max-parts",
                            "rna-reads inspect-concatemers",
                        )?;
                        settings.fragment_max_parts = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --fragment-max-parts value '{raw}' for rna-reads inspect-concatemers: {e}"
                            )
                        })?;
                    }
                    "--fragment-min-identity" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--fragment-min-identity",
                            "rna-reads inspect-concatemers",
                        )?;
                        settings.fragment_min_identity_fraction =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --fragment-min-identity value '{raw}' for rna-reads inspect-concatemers: {e}"
                                )
                            })?;
                    }
                    "--fragment-min-query-cov" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--fragment-min-query-cov",
                            "rna-reads inspect-concatemers",
                        )?;
                        settings.fragment_min_query_coverage_fraction =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --fragment-min-query-cov value '{raw}' for rna-reads inspect-concatemers: {e}"
                                )
                            })?;
                    }
                    "--transcript-fasta" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--transcript-fasta",
                            "rna-reads inspect-concatemers",
                        )?;
                        settings.transcript_fasta_paths.push(raw);
                    }
                    "--transcript-index" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--transcript-index",
                            "rna-reads inspect-concatemers",
                        )?;
                        settings.transcript_index_paths.push(raw);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads inspect-concatemers"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsInspectConcatemers {
                report_id,
                selection,
                limit,
                selected_record_indices,
                settings,
            })
        }
        "build-transcript-index" => {
            if tokens.len() < 4 {
                return Err(
                    "rna-reads build-transcript-index requires OUTPUT.json --transcript-fasta PATH [--transcript-fasta PATH ...] [--kmer-len N]"
                        .to_string(),
                );
            }
            let path = tokens[2].trim().to_string();
            if path.is_empty() {
                return Err(
                    "rna-reads build-transcript-index OUTPUT.json must not be empty".to_string(),
                );
            }
            let mut seed_kmer_len = 10usize;
            let mut transcript_fasta_paths = Vec::<String>::new();
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--transcript-fasta" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--transcript-fasta",
                            "rna-reads build-transcript-index",
                        )?;
                        transcript_fasta_paths.push(raw);
                    }
                    "--kmer-len" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--kmer-len",
                            "rna-reads build-transcript-index",
                        )?;
                        seed_kmer_len = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --kmer-len value '{raw}' for rna-reads build-transcript-index: {e}"
                            )
                        })?;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads build-transcript-index"
                        ));
                    }
                }
            }
            if transcript_fasta_paths.is_empty() {
                return Err(
                    "rna-reads build-transcript-index requires at least one --transcript-fasta PATH"
                        .to_string(),
                );
            }
            Ok(ShellCommand::RnaReadsBuildTranscriptIndex {
                path,
                seed_kmer_len,
                transcript_fasta_paths,
            })
        }
        "export-report" => {
            if tokens.len() != 4 {
                return Err("rna-reads export-report requires REPORT_ID OUTPUT.json".to_string());
            }
            Ok(ShellCommand::RnaReadsExportReport {
                report_id: tokens[2].clone(),
                path: tokens[3].clone(),
            })
        }
        "export-hits-fasta" => {
            if tokens.len() < 4 {
                return Err(
                    "rna-reads export-hits-fasta requires REPORT_ID OUTPUT.fa [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut selection = RnaReadHitSelection::Aligned;
            let mut selected_record_indices: Vec<usize> = vec![];
            let mut subset_spec: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--selection" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--selection",
                            "rna-reads export-hits-fasta",
                        )?;
                        selection = parse_rna_read_hit_selection(&raw)?;
                    }
                    "--record-indices" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--record-indices",
                            "rna-reads export-hits-fasta",
                        )?;
                        selected_record_indices = parse_rna_read_record_indices(&raw)?;
                    }
                    "--subset-spec" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--subset-spec",
                            "rna-reads export-hits-fasta",
                        )?;
                        subset_spec = Some(raw);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads export-hits-fasta"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsExportHitsFasta {
                report_id,
                path,
                selection,
                selected_record_indices,
                subset_spec,
            })
        }
        "export-sample-sheet" => {
            if tokens.len() < 3 {
                return Err(
                    "rna-reads export-sample-sheet requires OUTPUT.tsv [--seq-id ID] [--report-id ID]... [--gene GENE_ID]... [--complete-rule near|strict|exact] [--append]"
                        .to_string(),
                );
            }
            let path = tokens[2].clone();
            let mut seq_id: Option<String> = None;
            let mut report_ids: Vec<String> = vec![];
            let mut gene_ids: Vec<String> = vec![];
            let mut complete_rule = RnaReadGeneSupportCompleteRule::Near;
            let mut append = false;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--seq-id" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--seq-id",
                            "rna-reads export-sample-sheet",
                        )?;
                        seq_id = Some(raw);
                    }
                    "--report-id" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--report-id",
                            "rna-reads export-sample-sheet",
                        )?;
                        report_ids.push(raw);
                    }
                    "--gene" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--gene",
                            "rna-reads export-sample-sheet",
                        )?;
                        gene_ids.push(raw);
                    }
                    "--complete-rule" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--complete-rule",
                            "rna-reads export-sample-sheet",
                        )?;
                        complete_rule = parse_rna_read_gene_support_complete_rule(&raw)?;
                    }
                    "--append" => {
                        append = true;
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads export-sample-sheet"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsExportSampleSheet {
                path,
                seq_id,
                report_ids,
                gene_ids,
                complete_rule,
                append,
            })
        }
        "export-paths-tsv" => {
            if tokens.len() < 4 {
                return Err(
                    "rna-reads export-paths-tsv requires REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut selection = RnaReadHitSelection::All;
            let mut selected_record_indices: Vec<usize> = vec![];
            let mut subset_spec: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--selection" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--selection",
                            "rna-reads export-paths-tsv",
                        )?;
                        selection = parse_rna_read_hit_selection(&raw)?;
                    }
                    "--record-indices" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--record-indices",
                            "rna-reads export-paths-tsv",
                        )?;
                        selected_record_indices = parse_rna_read_record_indices(&raw)?;
                    }
                    "--subset-spec" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--subset-spec",
                            "rna-reads export-paths-tsv",
                        )?;
                        subset_spec = Some(raw);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads export-paths-tsv"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsExportExonPathsTsv {
                report_id,
                path,
                selection,
                selected_record_indices,
                subset_spec,
            })
        }
        "export-abundance-tsv" => {
            if tokens.len() < 4 {
                return Err(
                    "rna-reads export-abundance-tsv requires REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut selection = RnaReadHitSelection::All;
            let mut selected_record_indices: Vec<usize> = vec![];
            let mut subset_spec: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--selection" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--selection",
                            "rna-reads export-abundance-tsv",
                        )?;
                        selection = parse_rna_read_hit_selection(&raw)?;
                    }
                    "--record-indices" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--record-indices",
                            "rna-reads export-abundance-tsv",
                        )?;
                        selected_record_indices = parse_rna_read_record_indices(&raw)?;
                    }
                    "--subset-spec" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--subset-spec",
                            "rna-reads export-abundance-tsv",
                        )?;
                        subset_spec = Some(raw);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads export-abundance-tsv"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsExportExonAbundanceTsv {
                report_id,
                path,
                selection,
                selected_record_indices,
                subset_spec,
            })
        }
        "export-score-density-svg" => {
            if tokens.len() < 4 {
                return Err(
                    "rna-reads export-score-density-svg requires REPORT_ID OUTPUT.svg [--scale linear|log] [--variant all_scored|composite_seed_gate]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut scale = RnaReadScoreDensityScale::Log;
            let mut variant = RnaReadScoreDensityVariant::AllScored;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--scale" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--scale",
                            "rna-reads export-score-density-svg",
                        )?;
                        scale = parse_rna_read_score_density_scale(&raw)?;
                    }
                    "--variant" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--variant",
                            "rna-reads export-score-density-svg",
                        )?;
                        variant = parse_rna_read_score_density_variant(&raw)?;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads export-score-density-svg"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsExportScoreDensitySvg {
                report_id,
                path,
                scale,
                variant,
            })
        }
        "export-alignments-tsv" => {
            if tokens.len() < 4 {
                return Err(
                    "rna-reads export-alignments-tsv requires REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--limit N] [--record-indices i,j,k] [--subset-spec TEXT]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut selection = RnaReadHitSelection::Aligned;
            let mut limit: Option<usize> = None;
            let mut selected_record_indices: Vec<usize> = vec![];
            let mut subset_spec: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--selection" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--selection",
                            "rna-reads export-alignments-tsv",
                        )?;
                        selection = parse_rna_read_hit_selection(&raw)?;
                    }
                    "--limit" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--limit",
                            "rna-reads export-alignments-tsv",
                        )?;
                        limit = Some(raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --limit value '{raw}' for rna-reads export-alignments-tsv: {e}"
                            )
                        })?);
                    }
                    "--record-indices" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--record-indices",
                            "rna-reads export-alignments-tsv",
                        )?;
                        selected_record_indices = parse_rna_read_record_indices(&raw)?;
                    }
                    "--subset-spec" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--subset-spec",
                            "rna-reads export-alignments-tsv",
                        )?;
                        subset_spec = Some(raw);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads export-alignments-tsv"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsExportAlignmentsTsv {
                report_id,
                path,
                selection,
                limit,
                selected_record_indices,
                subset_spec,
            })
        }
        "export-alignment-dotplot-svg" => {
            if tokens.len() < 4 {
                return Err(
                    "rna-reads export-alignment-dotplot-svg requires REPORT_ID OUTPUT.svg [--selection all|seed_passed|aligned] [--max-points N]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut selection = RnaReadHitSelection::Aligned;
            let mut max_points = 2_500usize;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--selection" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--selection",
                            "rna-reads export-alignment-dotplot-svg",
                        )?;
                        selection = parse_rna_read_hit_selection(&raw)?;
                    }
                    "--max-points" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-points",
                            "rna-reads export-alignment-dotplot-svg",
                        )?;
                        max_points = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --max-points value '{raw}' for rna-reads export-alignment-dotplot-svg: {e}"
                            )
                        })?;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads export-alignment-dotplot-svg"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsExportAlignmentDotplotSvg {
                report_id,
                path,
                selection,
                max_points,
            })
        }
        other => Err(format!(
            "Unknown rna-reads subcommand '{other}' (expected interpret, align-report, list-reports, show-report, summarize-gene-support, inspect-gene-support, inspect-alignments, inspect-concatemers, build-transcript-index, export-report, export-hits-fasta, export-sample-sheet, export-paths-tsv, export-abundance-tsv, export-score-density-svg, export-alignments-tsv, export-alignment-dotplot-svg)"
        )),
    }
}

pub(super) fn parse_macros_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "macros requires a subcommand: run, instance-list, instance-show, template-list, template-show, template-put, template-delete, template-import, template-run"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "run" => {
            if tokens.len() < 3 {
                return Err(
                    "macros run requires SCRIPT_OR_@FILE (or --file PATH), optionally with --transactional".to_string(),
                );
            }
            let mut idx = 2usize;
            let mut transactional = false;
            let mut script_file: Option<String> = None;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--transactional" | "--atomic" => {
                        transactional = true;
                        idx += 1;
                    }
                    "--file" => {
                        if script_file.is_some() {
                            return Err("macros run --file may only be specified once".to_string());
                        }
                        idx += 1;
                        if idx >= tokens.len() {
                            return Err("macros run --file requires PATH".to_string());
                        }
                        script_file = Some(tokens[idx].trim().to_string());
                        idx += 1;
                    }
                    _ => break,
                }
            }
            let script = if let Some(path) = script_file {
                if idx != tokens.len() {
                    return Err(
                        "macros run does not accept inline script after --file PATH".to_string()
                    );
                }
                format!("@{path}")
            } else {
                if idx >= tokens.len() {
                    return Err("macros run requires SCRIPT_OR_@FILE".to_string());
                }
                tokens[idx..].join(" ")
            };
            if script.trim().is_empty() {
                return Err("macros run requires non-empty script".to_string());
            }
            Ok(ShellCommand::MacrosRun {
                script,
                transactional,
            })
        }
        "instance-list" => {
            if tokens.len() != 2 {
                return Err("macros instance-list takes no options".to_string());
            }
            Ok(ShellCommand::MacrosInstanceList)
        }
        "instance-show" => {
            if tokens.len() != 3 {
                return Err("macros instance-show requires MACRO_INSTANCE_ID".to_string());
            }
            Ok(ShellCommand::MacrosInstanceShow {
                macro_instance_id: tokens[2].clone(),
            })
        }
        "template-list" => {
            if tokens.len() != 2 {
                return Err("macros template-list takes no options".to_string());
            }
            Ok(ShellCommand::MacrosTemplateList)
        }
        "template-show" => {
            if tokens.len() != 3 {
                return Err("macros template-show requires TEMPLATE_NAME".to_string());
            }
            Ok(ShellCommand::MacrosTemplateShow {
                name: tokens[2].clone(),
            })
        }
        "template-put" | "template-upsert" => {
            if tokens.len() < 4 {
                return Err(
                    "macros template-put requires TEMPLATE_NAME (--script SCRIPT_OR_@FILE | --file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...] [--input-port PORT_ID:KIND[:one|many][:required|optional][:description]] [--output-port ...]"
                        .to_string(),
                );
            }
            let name = tokens[2].clone();
            let mut description: Option<String> = None;
            let mut details_url: Option<String> = None;
            let mut parameters: Vec<WorkflowMacroTemplateParam> = vec![];
            let mut input_ports: Vec<WorkflowMacroTemplatePort> = vec![];
            let mut output_ports: Vec<WorkflowMacroTemplatePort> = vec![];
            let mut script: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--description" => {
                        description = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--description",
                            "macros template-put",
                        )?);
                    }
                    "--details-url" | "--url" => {
                        if details_url.is_some() {
                            return Err(
                                "macros template-put details URL was already specified".to_string()
                            );
                        }
                        details_url = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--details-url",
                            "macros template-put",
                        )?);
                    }
                    "--param" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--param", "macros template-put")?;
                        parameters.push(parse_workflow_template_param_spec(&raw)?);
                    }
                    "--input-port" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--input-port",
                            "macros template-put",
                        )?;
                        input_ports.push(parse_workflow_template_port_spec(&raw)?);
                    }
                    "--output-port" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--output-port",
                            "macros template-put",
                        )?;
                        output_ports.push(parse_workflow_template_port_spec(&raw)?);
                    }
                    "--script" => {
                        if script.is_some() {
                            return Err(
                                "macros template-put script was already specified".to_string()
                            );
                        }
                        script = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--script",
                            "macros template-put",
                        )?);
                    }
                    "--file" => {
                        if script.is_some() {
                            return Err(
                                "macros template-put script was already specified".to_string()
                            );
                        }
                        let path =
                            parse_option_path(tokens, &mut idx, "--file", "macros template-put")?;
                        script = Some(format!("@{path}"));
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for macros template-put"));
                    }
                }
            }
            let script = script.ok_or_else(|| {
                "macros template-put requires --script SCRIPT_OR_@FILE or --file PATH".to_string()
            })?;
            Ok(ShellCommand::MacrosTemplateUpsert {
                name,
                description,
                details_url,
                parameters,
                input_ports,
                output_ports,
                script,
            })
        }
        "template-delete" => {
            if tokens.len() != 3 {
                return Err("macros template-delete requires TEMPLATE_NAME".to_string());
            }
            Ok(ShellCommand::MacrosTemplateDelete {
                name: tokens[2].clone(),
            })
        }
        "template-import" => {
            if tokens.len() != 3 {
                return Err("macros template-import requires PATH".to_string());
            }
            Ok(ShellCommand::MacrosTemplateImport {
                path: tokens[2].clone(),
            })
        }
        "template-run" => {
            if tokens.len() < 3 {
                return Err(
                    "macros template-run requires TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional] [--validate-only]"
                        .to_string(),
                );
            }
            let name = tokens[2].clone();
            let mut bindings: HashMap<String, String> = HashMap::new();
            let mut transactional = false;
            let mut validate_only = false;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--bind" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--bind", "macros template-run")?;
                        let (key, value) = parse_template_binding(&raw)?;
                        if bindings.insert(key.clone(), value).is_some() {
                            return Err(format!(
                                "Duplicate --bind key '{}' in macros template-run",
                                key
                            ));
                        }
                    }
                    "--transactional" | "--atomic" => {
                        transactional = true;
                        idx += 1;
                    }
                    "--validate-only" | "--dry-run" => {
                        validate_only = true;
                        idx += 1;
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for macros template-run"));
                    }
                }
            }
            Ok(ShellCommand::MacrosTemplateRun {
                name,
                bindings,
                transactional,
                validate_only,
            })
        }
        other => Err(format!(
            "Unknown macros subcommand '{other}' (expected run, instance-list, instance-show, template-list, template-show, template-put, template-delete, template-import, template-run)"
        )),
    }
}

pub(super) fn parse_routines_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("routines requires a subcommand: list, explain, compare".to_string());
    }
    match tokens[1].as_str() {
        "list" => {
            let mut catalog_path: Option<String> = None;
            let mut family: Option<String> = None;
            let mut status: Option<String> = None;
            let mut tag: Option<String> = None;
            let mut query: Option<String> = None;
            let mut seq_id: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--catalog",
                            "routines list",
                        )?);
                    }
                    "--family" => {
                        family = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--family",
                            "routines list",
                        )?);
                    }
                    "--status" => {
                        status = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--status",
                            "routines list",
                        )?);
                    }
                    "--tag" => {
                        tag = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--tag",
                            "routines list",
                        )?);
                    }
                    "--query" => {
                        query = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--query",
                            "routines list",
                        )?);
                    }
                    "--seq-id" => {
                        seq_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--seq-id",
                            "routines list",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for routines list"));
                    }
                }
            }
            Ok(ShellCommand::RoutinesList {
                catalog_path,
                family,
                status,
                tag,
                query,
                seq_id,
            })
        }
        "explain" => {
            if tokens.len() < 3 {
                return Err(
                    "routines explain requires ROUTINE_ID [--catalog PATH] [--seq-id SEQ_ID]"
                        .to_string(),
                );
            }
            let routine_id = tokens[2].trim().to_string();
            if routine_id.is_empty() {
                return Err("routines explain ROUTINE_ID cannot be empty".to_string());
            }
            let mut catalog_path: Option<String> = None;
            let mut seq_id: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--catalog",
                            "routines explain",
                        )?);
                    }
                    "--seq-id" => {
                        seq_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--seq-id",
                            "routines explain",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for routines explain"));
                    }
                }
            }
            Ok(ShellCommand::RoutinesExplain {
                catalog_path,
                routine_id,
                seq_id,
            })
        }
        "compare" => {
            if tokens.len() < 4 {
                return Err(
                    "routines compare requires ROUTINE_A ROUTINE_B [--catalog PATH] [--seq-id SEQ_ID]".to_string(),
                );
            }
            let left_routine_id = tokens[2].trim().to_string();
            let right_routine_id = tokens[3].trim().to_string();
            if left_routine_id.is_empty() || right_routine_id.is_empty() {
                return Err("routines compare routine ids cannot be empty".to_string());
            }
            let mut catalog_path: Option<String> = None;
            let mut seq_id: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--catalog",
                            "routines compare",
                        )?);
                    }
                    "--seq-id" => {
                        seq_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--seq-id",
                            "routines compare",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for routines compare"));
                    }
                }
            }
            Ok(ShellCommand::RoutinesCompare {
                catalog_path,
                left_routine_id,
                right_routine_id,
                seq_id,
            })
        }
        other => Err(format!(
            "Unknown routines subcommand '{other}' (expected list, explain, compare)"
        )),
    }
}

pub(super) fn parse_planning_profile_scope(raw: &str) -> Result<PlanningProfileScope, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "global" => Ok(PlanningProfileScope::Global),
        "project" | "project_override" | "project-override" => {
            Ok(PlanningProfileScope::ProjectOverride)
        }
        "agent" | "agent_overlay" | "confirmed_agent_overlay" | "confirmed-agent-overlay" => {
            Ok(PlanningProfileScope::ConfirmedAgentOverlay)
        }
        "effective" => Ok(PlanningProfileScope::Effective),
        other => Err(format!(
            "Unsupported planning profile scope '{other}' (expected global|project_override|confirmed_agent_overlay|effective)"
        )),
    }
}

pub(super) fn parse_planning_suggestion_status(
    raw: &str,
) -> Result<PlanningSuggestionStatus, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "pending" => Ok(PlanningSuggestionStatus::Pending),
        "accepted" => Ok(PlanningSuggestionStatus::Accepted),
        "rejected" => Ok(PlanningSuggestionStatus::Rejected),
        other => Err(format!(
            "Unsupported planning suggestion status '{other}' (expected pending|accepted|rejected)"
        )),
    }
}

pub(super) fn parse_planning_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "planning requires a subcommand: profile, objective, suggestions, sync".to_string(),
        );
    }
    match tokens[1].as_str() {
        "profile" => {
            if tokens.len() < 3 {
                return Err("planning profile requires a subcommand: show, set, clear".to_string());
            }
            match tokens[2].as_str() {
                "show" => {
                    let mut scope = PlanningProfileScope::Effective;
                    let mut idx = 3usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--scope" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--scope",
                                    "planning profile show",
                                )?;
                                scope = parse_planning_profile_scope(&raw)?;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for planning profile show"
                                ));
                            }
                        }
                    }
                    Ok(ShellCommand::PlanningProfileShow { scope })
                }
                "set" => {
                    if tokens.len() < 4 {
                        return Err(
                            "planning profile set requires JSON_OR_@FILE [--scope SCOPE]"
                                .to_string(),
                        );
                    }
                    let mut scope = PlanningProfileScope::ProjectOverride;
                    let payload_json = tokens[3].clone();
                    let mut idx = 4usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--scope" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--scope",
                                    "planning profile set",
                                )?;
                                scope = parse_planning_profile_scope(&raw)?;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for planning profile set"
                                ));
                            }
                        }
                    }
                    if scope == PlanningProfileScope::Effective {
                        return Err(
                            "planning profile set does not support --scope effective".to_string()
                        );
                    }
                    Ok(ShellCommand::PlanningProfileSet {
                        scope,
                        payload_json,
                    })
                }
                "clear" => {
                    let mut scope = PlanningProfileScope::ProjectOverride;
                    let mut idx = 3usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--scope" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--scope",
                                    "planning profile clear",
                                )?;
                                scope = parse_planning_profile_scope(&raw)?;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for planning profile clear"
                                ));
                            }
                        }
                    }
                    if scope == PlanningProfileScope::Effective {
                        return Err(
                            "planning profile clear does not support --scope effective".to_string()
                        );
                    }
                    Ok(ShellCommand::PlanningProfileSet {
                        scope,
                        payload_json: "null".to_string(),
                    })
                }
                other => Err(format!(
                    "Unknown planning profile subcommand '{other}' (expected show, set, clear)"
                )),
            }
        }
        "objective" => {
            if tokens.len() < 3 {
                return Err(
                    "planning objective requires a subcommand: show, set, clear".to_string()
                );
            }
            match tokens[2].as_str() {
                "show" => {
                    if tokens.len() > 3 {
                        return Err("planning objective show takes no options".to_string());
                    }
                    Ok(ShellCommand::PlanningObjectiveShow)
                }
                "set" => {
                    if tokens.len() != 4 {
                        return Err("planning objective set requires JSON_OR_@FILE".to_string());
                    }
                    Ok(ShellCommand::PlanningObjectiveSet {
                        payload_json: tokens[3].clone(),
                    })
                }
                "clear" => {
                    if tokens.len() > 3 {
                        return Err("planning objective clear takes no options".to_string());
                    }
                    Ok(ShellCommand::PlanningObjectiveSet {
                        payload_json: "null".to_string(),
                    })
                }
                other => Err(format!(
                    "Unknown planning objective subcommand '{other}' (expected show, set, clear)"
                )),
            }
        }
        "suggestions" => {
            if tokens.len() < 3 {
                return Err(
                    "planning suggestions requires a subcommand: list, accept, reject".to_string(),
                );
            }
            match tokens[2].as_str() {
                "list" => {
                    let mut status: Option<PlanningSuggestionStatus> = None;
                    let mut idx = 3usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--status" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--status",
                                    "planning suggestions list",
                                )?;
                                status = Some(parse_planning_suggestion_status(&raw)?);
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for planning suggestions list"
                                ));
                            }
                        }
                    }
                    Ok(ShellCommand::PlanningSuggestionsList { status })
                }
                "accept" => {
                    if tokens.len() != 4 {
                        return Err(
                            "planning suggestions accept requires SUGGESTION_ID".to_string()
                        );
                    }
                    let suggestion_id = tokens[3].trim().to_string();
                    if suggestion_id.is_empty() {
                        return Err(
                            "planning suggestions accept SUGGESTION_ID cannot be empty".to_string()
                        );
                    }
                    Ok(ShellCommand::PlanningSuggestionAccept { suggestion_id })
                }
                "reject" => {
                    if tokens.len() < 4 {
                        return Err(
                            "planning suggestions reject requires SUGGESTION_ID [--reason TEXT]"
                                .to_string(),
                        );
                    }
                    let suggestion_id = tokens[3].trim().to_string();
                    if suggestion_id.is_empty() {
                        return Err(
                            "planning suggestions reject SUGGESTION_ID cannot be empty".to_string()
                        );
                    }
                    let mut reason: Option<String> = None;
                    let mut idx = 4usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--reason" => {
                                reason = Some(parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--reason",
                                    "planning suggestions reject",
                                )?);
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for planning suggestions reject"
                                ));
                            }
                        }
                    }
                    Ok(ShellCommand::PlanningSuggestionReject {
                        suggestion_id,
                        reason,
                    })
                }
                other => Err(format!(
                    "Unknown planning suggestions subcommand '{other}' (expected list, accept, reject)"
                )),
            }
        }
        "sync" => {
            if tokens.len() < 3 {
                return Err("planning sync requires a subcommand: status, pull, push".to_string());
            }
            match tokens[2].as_str() {
                "status" => {
                    if tokens.len() > 3 {
                        return Err("planning sync status takes no options".to_string());
                    }
                    Ok(ShellCommand::PlanningSyncStatus)
                }
                "pull" | "push" => {
                    if tokens.len() < 4 {
                        return Err(format!(
                            "planning sync {} requires JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]",
                            tokens[2]
                        ));
                    }
                    let payload_json = tokens[3].clone();
                    let mut source: Option<String> = None;
                    let mut confidence: Option<f64> = None;
                    let mut snapshot_id: Option<String> = None;
                    let mut idx = 4usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--source" => {
                                source = Some(parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--source",
                                    "planning sync",
                                )?);
                            }
                            "--confidence" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--confidence",
                                    "planning sync",
                                )?;
                                let parsed = raw.parse::<f64>().map_err(|e| {
                                    format!("Invalid --confidence value '{raw}': {e}")
                                })?;
                                confidence = Some(parsed);
                            }
                            "--snapshot-id" => {
                                snapshot_id = Some(parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--snapshot-id",
                                    "planning sync",
                                )?);
                            }
                            other => {
                                return Err(format!("Unknown option '{other}' for planning sync"));
                            }
                        }
                    }
                    if tokens[2].eq_ignore_ascii_case("pull") {
                        Ok(ShellCommand::PlanningSyncPull {
                            payload_json,
                            source,
                            confidence,
                            snapshot_id,
                        })
                    } else {
                        Ok(ShellCommand::PlanningSyncPush {
                            payload_json,
                            source,
                            confidence,
                            snapshot_id,
                        })
                    }
                }
                other => Err(format!(
                    "Unknown planning sync subcommand '{other}' (expected status, pull, push)"
                )),
            }
        }
        other => Err(format!(
            "Unknown planning subcommand '{other}' (expected profile, objective, suggestions, sync)"
        )),
    }
}
