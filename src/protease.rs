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
