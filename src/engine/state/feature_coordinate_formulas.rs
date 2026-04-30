//! Shared feature-relative coordinate formula parsing helpers.
//!
//! This module centralizes the `=KIND.start+N` / `=left .. right` parsing used
//! by GUI selection and ROI-entry controls so future shell/CLI surfaces can
//! reuse the exact same resolution logic instead of re-implementing it in
//! adapters.

use crate::{dna_sequence::DNAsequence, feature_location::collect_location_ranges_usize};

use super::AnchorBoundary;

fn feature_formula_label_values(feature: &gb_io::seq::Feature) -> Vec<String> {
    let mut labels = Vec::new();
    for key in [
        "label",
        "gene",
        "locus_tag",
        "product",
        "standard_name",
        "note",
    ] {
        for value in feature.qualifier_values(key.into()) {
            let value = value.trim();
            if !value.is_empty() {
                labels.push(value.to_string());
            }
        }
    }
    labels
}

pub fn split_feature_formula_range_expression(expr: &str) -> Option<(String, String)> {
    if let Some((left, right)) = expr.split_once("..") {
        let left = left.trim();
        let right = right.trim();
        if !left.is_empty() && !right.is_empty() {
            return Some((left.to_string(), right.to_string()));
        }
    }
    let lower = expr.to_ascii_lowercase();
    if let Some(split_idx) = lower.find(" to ") {
        let left = expr[..split_idx].trim();
        let right = expr[split_idx + 4..].trim();
        if !left.is_empty() && !right.is_empty() {
            return Some((left.to_string(), right.to_string()));
        }
    }
    None
}

fn parse_feature_formula_coordinate_expression_on_sequence(
    dna: &DNAsequence,
    expr: &str,
    field_name: &str,
) -> Result<usize, String> {
    let raw = expr.trim();
    if raw.is_empty() {
        return Err(format!(
            "Invalid {field_name}: empty feature formula expression"
        ));
    }

    let mut idx = 0usize;
    let bytes = raw.as_bytes();
    let is_kind_char = |ch: u8| ch.is_ascii_alphanumeric() || matches!(ch, b'_' | b'-');
    while idx < bytes.len() && is_kind_char(bytes[idx]) {
        idx += 1;
    }
    if idx == 0 {
        return Err(format!(
            "Invalid {field_name}: expected feature kind (for example CDS.start+10)"
        ));
    }
    let feature_kind = raw[..idx].trim().to_ascii_uppercase();

    while idx < bytes.len() && bytes[idx].is_ascii_whitespace() {
        idx += 1;
    }

    let mut occurrence_index: Option<usize> = None;
    let mut label_filter: Option<String> = None;
    if idx < bytes.len() && bytes[idx] == b'[' {
        idx += 1;
        let start = idx;
        while idx < bytes.len() && bytes[idx] != b']' {
            idx += 1;
        }
        if idx >= bytes.len() || bytes[idx] != b']' {
            return Err(format!(
                "Invalid {field_name}: missing closing ']' in feature selector"
            ));
        }
        let inside = raw[start..idx].trim();
        idx += 1;
        if inside.is_empty() {
            return Err(format!("Invalid {field_name}: empty selector inside []"));
        }
        let inside_lower = inside.to_ascii_lowercase();
        if let Some(value_raw) = inside_lower.strip_prefix("label=").map(|_| &inside[6..]) {
            let label = value_raw.trim();
            if label.is_empty() {
                return Err(format!(
                    "Invalid {field_name}: label filter must not be empty"
                ));
            }
            label_filter = Some(label.to_ascii_uppercase());
        } else {
            let one_based = inside.parse::<usize>().map_err(|_| {
                format!(
                    "Invalid {field_name}: selector '{}' must be a 1-based occurrence index or label=...",
                    inside
                )
            })?;
            if one_based == 0 {
                return Err(format!(
                    "Invalid {field_name}: occurrence index must be >= 1"
                ));
            }
            occurrence_index = Some(one_based - 1);
        }
    }

    let mut saw_separator = false;
    if idx < bytes.len() && bytes[idx] == b'.' {
        saw_separator = true;
        idx += 1;
    }
    while idx < bytes.len() && bytes[idx].is_ascii_whitespace() {
        saw_separator = true;
        idx += 1;
    }
    if !saw_separator {
        return Err(format!(
            "Invalid {field_name}: expected `.start`, `.end`, or `.middle` after feature kind"
        ));
    }

    let boundary_start = idx;
    while idx < bytes.len() && bytes[idx].is_ascii_alphabetic() {
        idx += 1;
    }
    if boundary_start == idx {
        return Err(format!(
            "Invalid {field_name}: expected boundary token (start|end|middle)"
        ));
    }
    let boundary = match raw[boundary_start..idx]
        .trim()
        .to_ascii_lowercase()
        .as_str()
    {
        "start" => AnchorBoundary::Start,
        "end" => AnchorBoundary::End,
        "middle" => AnchorBoundary::Middle,
        other => {
            return Err(format!(
                "Invalid {field_name}: unknown boundary '{other}' (expected start|end|middle)"
            ));
        }
    };

    let mut offset: isize = 0;
    while idx < bytes.len() {
        while idx < bytes.len() && bytes[idx].is_ascii_whitespace() {
            idx += 1;
        }
        if idx >= bytes.len() {
            break;
        }
        let sign = match bytes[idx] {
            b'+' => 1isize,
            b'-' => -1isize,
            other => {
                return Err(format!(
                    "Invalid {field_name}: unexpected character '{}' in offset suffix",
                    other as char
                ));
            }
        };
        idx += 1;
        while idx < bytes.len() && bytes[idx].is_ascii_whitespace() {
            idx += 1;
        }
        let number_start = idx;
        while idx < bytes.len() && bytes[idx].is_ascii_digit() {
            idx += 1;
        }
        if number_start == idx {
            return Err(format!(
                "Invalid {field_name}: expected integer offset after sign"
            ));
        }
        let delta = raw[number_start..idx]
            .parse::<isize>()
            .map_err(|_| format!("Invalid {field_name}: could not parse coordinate offset"))?;
        offset += sign * delta;
    }

    if dna.len() == 0 {
        return Err(format!("Invalid {field_name}: active sequence is empty"));
    }

    let mut matches: Vec<(usize, usize, usize)> = Vec::new();
    for (feature_id, feature) in dna.features().iter().enumerate() {
        if feature.kind.to_string().eq_ignore_ascii_case("SOURCE")
            || !feature.kind.to_string().eq_ignore_ascii_case(&feature_kind)
        {
            continue;
        }
        if let Some(label_filter) = label_filter.as_ref() {
            let labels = feature_formula_label_values(feature);
            let found = labels.iter().any(|label| {
                let upper = label.to_ascii_uppercase();
                upper == *label_filter || upper.contains(label_filter)
            });
            if !found {
                continue;
            }
        }
        let mut ranges: Vec<(usize, usize)> = Vec::new();
        collect_location_ranges_usize(&feature.location, &mut ranges);
        if ranges.is_empty() {
            let Ok((from, to)) = feature.location.find_bounds() else {
                continue;
            };
            if from < 0 || to < 0 {
                continue;
            }
            ranges.push((from as usize, to as usize));
        }
        let start = match ranges.iter().map(|(start, _)| *start).min() {
            Some(value) => value,
            None => continue,
        };
        let end_exclusive = match ranges.iter().map(|(_, end)| *end).max() {
            Some(value) => value.min(dna.len()),
            None => continue,
        };
        if end_exclusive <= start || start >= dna.len() {
            continue;
        }
        let anchor_pos = match boundary {
            AnchorBoundary::Start => start,
            AnchorBoundary::End => end_exclusive,
            AnchorBoundary::Middle => start + (end_exclusive.saturating_sub(start) / 2),
        };
        if anchor_pos <= dna.len() {
            matches.push((start, feature_id, anchor_pos));
        }
    }

    if matches.is_empty() {
        return Err(format!(
            "Invalid {field_name}: no feature matched {}{}",
            feature_kind,
            label_filter
                .as_ref()
                .map(|label| format!("[label={label}]"))
                .unwrap_or_default()
        ));
    }
    matches.sort_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));
    let occurrence_idx = occurrence_index.unwrap_or(0);
    let Some((_, _, base_pos)) = matches.get(occurrence_idx) else {
        return Err(format!(
            "Invalid {field_name}: occurrence {} requested but only {} match(es) found",
            occurrence_idx + 1,
            matches.len()
        ));
    };
    let resolved = *base_pos as isize + offset;
    if resolved < 0 || resolved > dna.len() as isize {
        return Err(format!(
            "Invalid {field_name}: resolved coordinate {} is out of bounds for sequence length {}",
            resolved,
            dna.len()
        ));
    }
    Ok(resolved as usize)
}

pub fn parse_feature_coordinate_term_on_sequence(
    dna: &DNAsequence,
    raw: &str,
    field_name: &str,
) -> Result<usize, String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err(format!(
            "Invalid {field_name}: expected an integer or formula"
        ));
    }
    if let Ok(value) = trimmed.parse::<usize>() {
        return Ok(value);
    }
    parse_feature_formula_coordinate_expression_on_sequence(dna, trimmed, field_name)
}

pub fn parse_required_usize_or_formula_text_on_sequence(
    dna: &DNAsequence,
    raw: &str,
    field_name: &str,
) -> Result<usize, String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err(format!("Invalid {field_name}: expected an integer"));
    }
    if let Some(expr) = trimmed.strip_prefix('=') {
        return parse_feature_coordinate_term_on_sequence(dna, expr, field_name);
    }
    trimmed
        .parse::<usize>()
        .map_err(|_| format!("Invalid {field_name}: expected an integer"))
}

pub fn resolve_formula_roi_range_inputs_0based_on_sequence(
    dna: &DNAsequence,
    roi_start_raw: &str,
    roi_end_raw: &str,
    field_prefix: &str,
) -> Result<(usize, usize), String> {
    let start_trimmed = roi_start_raw.trim();
    let (start, end_exclusive) = if let Some(formula) = start_trimmed.strip_prefix('=')
        && let Some((left, right)) = split_feature_formula_range_expression(formula)
    {
        (
            parse_feature_coordinate_term_on_sequence(
                dna,
                &left,
                &format!("{field_prefix}.roi_start_0based"),
            )?,
            parse_feature_coordinate_term_on_sequence(
                dna,
                &right,
                &format!("{field_prefix}.roi_end_0based"),
            )?,
        )
    } else {
        (
            parse_required_usize_or_formula_text_on_sequence(
                dna,
                roi_start_raw,
                &format!("{field_prefix}.roi_start_0based"),
            )?,
            parse_required_usize_or_formula_text_on_sequence(
                dna,
                roi_end_raw,
                &format!("{field_prefix}.roi_end_0based"),
            )?,
        )
    };
    if end_exclusive <= start {
        return Err(format!(
            "Invalid {field_prefix} ROI range: start ({start}) must be < end ({end_exclusive})"
        ));
    }
    let seq_len = dna.len();
    if seq_len == 0 {
        return Err(format!("Invalid {field_prefix}: active sequence is empty"));
    }
    if start >= seq_len {
        return Err(format!(
            "Invalid {field_prefix}.roi_start_0based: {start} is outside sequence length {seq_len}"
        ));
    }
    if end_exclusive > seq_len {
        return Err(format!(
            "Invalid {field_prefix}.roi_end_0based: {end_exclusive} is outside sequence length {seq_len}"
        ));
    }
    Ok((start, end_exclusive))
}

pub fn resolve_selection_formula_range_0based_on_sequence(
    dna: &DNAsequence,
    raw: &str,
) -> Result<(usize, usize), String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err(
            "Selection formula is empty; expected `=left .. right` or `=left to right`".to_string(),
        );
    }
    let expr = trimmed
        .strip_prefix('=')
        .ok_or_else(|| "Selection formula must start with '='".to_string())?;
    let (left_raw, right_raw) = split_feature_formula_range_expression(expr).ok_or_else(|| {
        "Selection formula must define a range (`=left .. right` or `=left to right`)".to_string()
    })?;
    let start =
        parse_feature_coordinate_term_on_sequence(dna, &left_raw, "selection_formula.start")?;
    let end_exclusive =
        parse_feature_coordinate_term_on_sequence(dna, &right_raw, "selection_formula.end")?;
    if end_exclusive <= start {
        return Err(format!(
            "Invalid selection formula range: start ({start}) must be < end ({end_exclusive})"
        ));
    }
    let seq_len = dna.len();
    if seq_len == 0 {
        return Err("Cannot apply selection formula: active sequence is empty".to_string());
    }
    if start >= seq_len {
        return Err(format!(
            "Invalid selection formula start: {start} is outside sequence length {seq_len}"
        ));
    }
    if end_exclusive > seq_len {
        return Err(format!(
            "Invalid selection formula end: {end_exclusive} is outside sequence length {seq_len}"
        ));
    }
    Ok((start, end_exclusive))
}
