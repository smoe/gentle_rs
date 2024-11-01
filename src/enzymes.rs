use std::fs;

use crate::{error::GENtleError, protease::Protease, restriction_enzyme::RestrictionEnzyme};

#[derive(Clone, Debug, Default)]
pub struct Enzymes {
    restriction_enzymes: Vec<RestrictionEnzyme>,
    proteases: Vec<Protease>,
    max_re_length: usize,
    has_nonpalindromic_restriction_enzymes: bool,
}

impl Enzymes {
    pub fn from_json_file(filename: &str) -> Result<Self, GENtleError> {
        let mut ret = Self::default();
        let data = fs::read_to_string(filename)?;
        let res: serde_json::Value = serde_json::from_str(&data)?;
        let arr = res.as_array().ok_or(GENtleError::String(format!(
            "JSON in file '{filename}' is not an array"
        )))?;
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
                            None => {
                                return Err(GENtleError::String(format!(
                                    "Bad restriction enzyme: {row}"
                                )))
                            }
                        };
                    re.check_palimdromic();
                    ret.restriction_enzymes.push(re);
                }
                Some("protease") => {
                    let pt: Protease = match serde_json::from_str(&row.to_string()).ok() {
                        Some(pt) => pt,
                        None => return Err(GENtleError::String(format!("Bad protease: {row}"))),
                    };
                    ret.proteases.push(pt);
                }
                Some(other) => {
                    return Err(GENtleError::String(format!(
                        "Unknown enzyme type '{other}' in {}",
                        row
                    )))
                }
                None => {
                    return Err(GENtleError::String(format!(
                        "Missing enzyme type for {}",
                        row
                    )))
                }
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
        let enzymes = Enzymes::from_json_file("assets/enzymes.json").unwrap();
        assert!(enzymes
            .restriction_enzymes
            .iter()
            .any(|e| e.name == "EcoRI"));
        assert!(enzymes.proteases.iter().any(|e| e.name == "Clostripain"));
    }
}
