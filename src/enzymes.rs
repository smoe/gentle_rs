//! Restriction-enzyme catalog loading and convenience selection helpers.

use crate::{protease::Protease, restriction_enzyme::RestrictionEnzyme};
use anyhow::{Result, anyhow};
use std::fs;

const RUNTIME_REBASE_PATH: &str = "data/resources/rebase.enzymes.json";
const BUILTIN_ENZYMES_JSON: &str = include_str!("../assets/enzymes.json");
const DEFAULT_PREFERRED_RESTRICTION_ENZYME_NAMES: &[&str] = &[
    "EcoRI", "SacI", "KpnI", "SmaI", "BamHI", "XbaI", "SalI", "PstI", "SphI", "HindIII",
];
const GOLDEN_GATE_TYPE_IIS_RESTRICTION_ENZYME_NAMES: &[&str] = &[
    "BsaI", "Eco31I", "Eco31", "BsmBI", "Esp3I", "BbsI", "AarI", "SapI", "BtgZI", "BsmAI", "BfuAI",
];

#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct Enzymes {
    restriction_enzymes: Vec<RestrictionEnzyme>,
    proteases: Vec<Protease>,
    max_re_length: usize,
    has_nonpalindromic_restriction_enzymes: bool,
}

impl Enzymes {
    fn new(json_text: &str) -> Result<Self> {
        let mut ret = Self {
            restriction_enzymes: vec![],
            proteases: vec![],
            max_re_length: 0,
            has_nonpalindromic_restriction_enzymes: false,
        };
        let res: serde_json::Value = serde_json::from_str(json_text)?;
        let arr = res
            .as_array()
            .ok_or(anyhow!("Enzymes file is not a JSON array"))?;
        for row in arr {
            let enzyme_type = match row.get("type") {
                Some(et) => et,
                None => continue, // TODO warning?
            };
            match enzyme_type.as_str() {
                Some("restriction") => {
                    let mut re: RestrictionEnzyme =
                        match serde_json::from_str(&row.to_string()).ok() {
                            Some(re) => re,
                            None => return Err(anyhow!("Bad restriction enzyme: {row}")),
                        };
                    re.check_palimdromic();
                    ret.restriction_enzymes.push(re);
                }
                Some("protease") => {
                    let pt: Protease = match serde_json::from_str(&row.to_string()).ok() {
                        Some(pt) => pt,
                        None => return Err(anyhow!("Bad protease: {row}")),
                    };
                    ret.proteases.push(pt);
                }
                Some(other) => return Err(anyhow!("Unknown enzyme type '{other}' in {}", row)),
                None => return Err(anyhow!("Missing enzyme type for {}", row)),
            }
        }
        ret.max_re_length = ret
            .restriction_enzymes
            .iter()
            .map(|re| re.sequence.len())
            .max()
            .unwrap_or(0);
        ret.has_nonpalindromic_restriction_enzymes = ret
            .restriction_enzymes
            .iter()
            .any(|re| !re.is_palindromic());
        Ok(ret)
    }

    pub fn restriction_enzymes(&self) -> &Vec<RestrictionEnzyme> {
        &self.restriction_enzymes
    }

    pub fn proteases(&self) -> &Vec<Protease> {
        &self.proteases
    }

    pub fn restriction_enzymes_by_name(&self, names: &[&str]) -> Vec<RestrictionEnzyme> {
        self.restriction_enzymes
            .iter()
            .filter(|re| names.contains(&re.name.as_str()))
            .cloned()
            .collect()
    }

    fn recompute_derived_fields(&mut self) {
        self.max_re_length = self
            .restriction_enzymes
            .iter()
            .map(|re| re.sequence.len())
            .max()
            .unwrap_or(0);
        self.has_nonpalindromic_restriction_enzymes = self
            .restriction_enzymes
            .iter()
            .any(|re| !re.is_palindromic());
    }
}

pub fn load_restriction_enzymes_from_json_text(json_text: &str) -> Result<Vec<RestrictionEnzyme>> {
    Ok(Enzymes::new(json_text)?.restriction_enzymes)
}

pub fn load_restriction_enzymes_from_path(path: &str) -> Result<Vec<RestrictionEnzyme>> {
    let text = fs::read_to_string(path)?;
    load_restriction_enzymes_from_json_text(&text)
}

pub fn active_restriction_enzymes() -> Vec<RestrictionEnzyme> {
    if let Ok(custom) = load_restriction_enzymes_from_path(RUNTIME_REBASE_PATH) {
        if !custom.is_empty() {
            return custom;
        }
    }
    load_restriction_enzymes_from_json_text(BUILTIN_ENZYMES_JSON).unwrap_or_default()
}

pub fn default_preferred_restriction_enzyme_names() -> Vec<String> {
    DEFAULT_PREFERRED_RESTRICTION_ENZYME_NAMES
        .iter()
        .map(|name| (*name).to_string())
        .collect()
}

pub fn normalize_restriction_enzyme_name_token(raw: &str) -> String {
    raw.chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .map(|c| c.to_ascii_lowercase())
        .collect()
}

pub fn is_type_iis_capable_enzyme_name(name: &str) -> bool {
    matches!(
        normalize_restriction_enzyme_name_token(name).as_str(),
        "eco31"
            | "eco31i"
            | "bsai"
            | "bsmbi"
            | "esp3i"
            | "bbsi"
            | "aari"
            | "sapi"
            | "btgzi"
            | "bsmai"
            | "bfuai"
    )
}

pub fn golden_gate_type_iis_preferred_restriction_enzyme_names_from_catalog(
    enzymes: &[RestrictionEnzyme],
) -> Vec<String> {
    let mut lookup = std::collections::BTreeMap::new();
    for enzyme in enzymes {
        lookup.insert(
            normalize_restriction_enzyme_name_token(&enzyme.name),
            enzyme.name.clone(),
        );
    }
    let mut selected = vec![];
    for candidate in GOLDEN_GATE_TYPE_IIS_RESTRICTION_ENZYME_NAMES {
        let normalized = normalize_restriction_enzyme_name_token(candidate);
        if let Some(actual_name) = lookup.remove(&normalized) {
            selected.push(actual_name);
        }
    }
    selected
}

pub fn golden_gate_type_iis_preferred_restriction_enzyme_names() -> Vec<String> {
    let enzymes = active_restriction_enzymes();
    golden_gate_type_iis_preferred_restriction_enzyme_names_from_catalog(&enzymes)
}

impl Default for Enzymes {
    fn default() -> Self {
        let mut base = Enzymes::new(BUILTIN_ENZYMES_JSON).unwrap();
        if let Ok(text) = fs::read_to_string(RUNTIME_REBASE_PATH) {
            if let Ok(custom) = Enzymes::new(&text) {
                if !custom.restriction_enzymes.is_empty() {
                    base.restriction_enzymes = custom.restriction_enzymes;
                    base.recompute_derived_fields();
                }
            }
        }
        base
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_json_file() {
        let enzymes = Enzymes::default();
        assert!(
            enzymes
                .restriction_enzymes
                .iter()
                .any(|e| e.name == "EcoRI")
        );
        assert!(enzymes.proteases.iter().any(|e| e.name == "Clostripain"));
    }

    #[test]
    fn golden_gate_type_iis_preset_prefers_common_catalog_names() {
        let mk = |name: &str, sequence: &str, cut: isize, overlap: isize| {
            serde_json::from_value::<RestrictionEnzyme>(serde_json::json!({
                "name": name,
                "sequence": sequence,
                "note": null,
                "cut": cut,
                "overlap": overlap
            }))
            .expect("synthetic restriction enzyme")
        };
        let enzymes = vec![
            mk("EcoRI", "GAATTC", 1, 0),
            mk("BsaI", "GGTCTC", 7, -4),
            mk("Esp3I", "CGTCTC", 7, -4),
            mk("BbsI", "GAAGAC", 8, -4),
        ];

        assert_eq!(
            golden_gate_type_iis_preferred_restriction_enzyme_names_from_catalog(&enzymes),
            vec!["BsaI".to_string(), "Esp3I".to_string(), "BbsI".to_string()]
        );
        assert!(is_type_iis_capable_enzyme_name("BsmBI"));
        assert!(is_type_iis_capable_enzyme_name("Eco31"));
        assert!(!is_type_iis_capable_enzyme_name("EcoRI"));
    }
}
