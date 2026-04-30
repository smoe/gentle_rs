//! Protease digest definitions and helpers.

use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct Protease {
    pub name: String,
    pub sequence: String,
    pub note: Option<String>,
    pub cut: isize,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub aliases: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub category: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub specificity: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cleavage_side: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub typical_applications: Vec<String>,
    #[serde(default)]
    pub sequence_specific: bool,
}

pub fn normalize_protease_name_token(raw: &str) -> String {
    raw.chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .map(|c| c.to_ascii_lowercase())
        .collect()
}

impl Protease {
    pub fn cleavage_boundaries_0based(&self, protein_sequence: &str) -> Vec<usize> {
        let Some((left, right)) = self.sequence.split_once('|') else {
            return vec![];
        };
        let left_tokens = split_protease_pattern_tokens(left);
        let right_tokens = split_protease_pattern_tokens(right);
        if left_tokens.is_empty() && right_tokens.is_empty() {
            return vec![];
        }
        let residues = protein_sequence
            .as_bytes()
            .iter()
            .map(|residue| residue.to_ascii_uppercase())
            .collect::<Vec<_>>();
        let len = residues.len();
        let mut boundaries = vec![];
        for boundary in 0..=len {
            if boundary < left_tokens.len() || boundary + right_tokens.len() > len {
                continue;
            }
            let left_start = boundary.saturating_sub(left_tokens.len());
            let left_matches = left_tokens.iter().enumerate().all(|(idx, token)| {
                protease_pattern_token_matches(token, residues[left_start + idx])
            });
            if !left_matches {
                continue;
            }
            let right_matches = right_tokens.iter().enumerate().all(|(idx, token)| {
                protease_pattern_token_matches(token, residues[boundary + idx])
            });
            if right_matches {
                boundaries.push(boundary);
            }
        }
        boundaries
    }

    pub fn matches_filter(&self, filter: &str) -> bool {
        let needle = filter.trim().to_ascii_lowercase();
        if needle.is_empty() {
            return true;
        }
        self.search_haystacks()
            .into_iter()
            .any(|value| value.to_ascii_lowercase().contains(&needle))
    }

    pub fn matches_name_or_alias(&self, query: &str) -> bool {
        let needle = normalize_protease_name_token(query);
        if needle.is_empty() {
            return false;
        }
        normalize_protease_name_token(&self.name) == needle
            || self
                .aliases
                .iter()
                .any(|alias| normalize_protease_name_token(alias) == needle)
    }

    pub fn search_haystacks(&self) -> Vec<&str> {
        let mut haystacks = vec![self.name.as_str(), self.sequence.as_str()];
        if let Some(note) = self.note.as_deref() {
            haystacks.push(note);
        }
        if let Some(category) = self.category.as_deref() {
            haystacks.push(category);
        }
        if let Some(specificity) = self.specificity.as_deref() {
            haystacks.push(specificity);
        }
        if let Some(cleavage_side) = self.cleavage_side.as_deref() {
            haystacks.push(cleavage_side);
        }
        for alias in &self.aliases {
            haystacks.push(alias.as_str());
        }
        for application in &self.typical_applications {
            haystacks.push(application.as_str());
        }
        haystacks
    }
}

fn split_protease_pattern_tokens(pattern: &str) -> Vec<String> {
    let trimmed = pattern.trim();
    if trimmed.is_empty() {
        return vec![];
    }
    if trimmed.contains(',') {
        trimmed
            .split(',')
            .map(str::trim)
            .filter(|token| !token.is_empty())
            .map(ToOwned::to_owned)
            .collect()
    } else {
        vec![trimmed.to_string()]
    }
}

fn protease_pattern_token_matches(token: &str, residue: u8) -> bool {
    let token = token.trim().to_ascii_uppercase();
    if token.is_empty() {
        return false;
    }
    if token == "*" {
        return true;
    }
    if let Some(excluded) = token.strip_prefix('!') {
        return !excluded.as_bytes().contains(&residue);
    }
    token.as_bytes().contains(&residue)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn protease(name: &str, pattern: &str) -> Protease {
        Protease {
            name: name.to_string(),
            sequence: pattern.to_string(),
            ..Protease::default()
        }
    }

    #[test]
    fn trypsin_boundaries_skip_proline_blocked_sites() {
        let trypsin = protease("Trypsin", "KR|!P");
        assert_eq!(
            trypsin.cleavage_boundaries_0based("MAKRPTRKAA"),
            vec![3, 7, 8]
        );
    }

    #[test]
    fn position_specific_tag_protease_patterns_match_boundary() {
        let tev = protease("TEV", "E,N,L,Y,F,Q|G");
        assert_eq!(
            tev.cleavage_boundaries_0based("AAENLYFQGSSENLYFQASS"),
            vec![8]
        );
    }

    #[test]
    fn n_terminal_specificity_patterns_match_before_target_residue() {
        let lys_n = protease("Lys-N", "*|K");
        assert_eq!(lys_n.cleavage_boundaries_0based("MAKTEKAA"), vec![2, 5]);
    }
}
