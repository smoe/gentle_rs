use crate::{protease::Protease, restriction_enzyme::RestrictionEnzyme};
use anyhow::{anyhow, Result};

#[derive(Clone, Debug, Default)]
pub struct Enzymes {
    restriction_enzymes: Vec<RestrictionEnzyme>,
    proteases: Vec<Protease>,
    max_re_length: usize,
    has_nonpalindromic_restriction_enzymes: bool,
}

impl Enzymes {
    pub fn new() -> Result<Self> {
        let mut ret = Self::default();
        let data = include_str!("../assets/enzymes.json");
        let res: serde_json::Value = serde_json::from_str(data)?;
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_json_file() {
        let enzymes = Enzymes::new().unwrap();
        assert!(enzymes
            .restriction_enzymes
            .iter()
            .any(|e| e.name == "EcoRI"));
        assert!(enzymes.proteases.iter().any(|e| e.name == "Clostripain"));
    }
}
