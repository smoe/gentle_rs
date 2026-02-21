use crate::{protease::Protease, restriction_enzyme::RestrictionEnzyme};
use anyhow::{Result, anyhow};
use std::fs;

const RUNTIME_REBASE_PATH: &str = "data/resources/rebase.enzymes.json";
const BUILTIN_ENZYMES_JSON: &str = include_str!("../assets/enzymes.json");

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
}
