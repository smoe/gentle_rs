use serde_json::Value;
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum LadderMolecule {
    Dna,
    Rna,
}

impl LadderMolecule {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Dna => "dna",
            Self::Rna => "rna",
        }
    }

    pub fn display_name(self) -> &'static str {
        match self {
            Self::Dna => "DNA",
            Self::Rna => "RNA",
        }
    }

    pub fn parse(text: &str) -> Option<Self> {
        let norm = text.trim().to_ascii_lowercase();
        match norm.as_str() {
            "dna" => Some(Self::Dna),
            "rna" => Some(Self::Rna),
            _ => None,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct LadderBand {
    length: f64,
    pub relative_strength: Option<f64>,
}

impl LadderBand {
    pub fn length(&self) -> f64 {
        self.length
    }

    pub fn length_bp(&self) -> f64 {
        self.length
    }

    pub fn length_nt(&self) -> f64 {
        self.length
    }
}

#[derive(Clone, Debug, Default)]
pub struct Ladder {
    name: String,
    loading_hint: Option<f64>,
    bands: Vec<LadderBand>,
}

impl Ladder {
    pub fn new(name: &str, parts: &Value) -> Self {
        let raw_parts = parts
            .as_array()
            .expect("DNA ladder part is not an array")
            .iter()
            .map(|p| p.as_array().expect("DNA ladder subpart not an array"))
            .collect::<Vec<_>>();

        // Historical ladder files often store a non-band load hint as first entry.
        let (start_idx, loading_hint) = if raw_parts.len() >= 2
            && raw_parts[0].len() == 1
            && raw_parts[0][0].is_number()
            && raw_parts[1].first().map(|v| v.is_number()).unwrap_or(false)
        {
            let first = raw_parts[0][0].as_f64().unwrap_or_default();
            let second = raw_parts[1][0].as_f64().unwrap_or_default();
            if first > 0.0 && second > first * 2.0 {
                (1, Some(first))
            } else {
                (0, None)
            }
        } else {
            (0, None)
        };

        let bands: Vec<LadderBand> = raw_parts
            .into_iter()
            .skip(start_idx)
            .filter_map(|p| {
                let length = p.first()?.as_f64()?;
                if !length.is_finite() || length <= 0.0 {
                    return None;
                }
                Some(LadderBand {
                    length,
                    relative_strength: p.get(1).and_then(|s| s.as_f64()),
                })
            })
            .collect();

        Self {
            name: name.to_owned(),
            loading_hint,
            bands,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn loading_hint(&self) -> Option<f64> {
        self.loading_hint
    }

    pub fn bands(&self) -> &Vec<LadderBand> {
        &self.bands
    }

    pub fn min_len(&self) -> Option<usize> {
        self.bands.iter().map(|p| p.length.round() as usize).min()
    }

    pub fn max_len(&self) -> Option<usize> {
        self.bands.iter().map(|p| p.length.round() as usize).max()
    }

    pub fn min_bp(&self) -> Option<usize> {
        self.min_len()
    }

    pub fn max_bp(&self) -> Option<usize> {
        self.max_len()
    }

    pub fn min_nt(&self) -> Option<usize> {
        self.min_len()
    }

    pub fn max_nt(&self) -> Option<usize> {
        self.max_len()
    }
}

#[derive(Clone, Debug)]
pub struct LadderCatalog {
    molecule: LadderMolecule,
    ladders: HashMap<String, Ladder>,
}

impl LadderCatalog {
    pub fn from_json_str(molecule: LadderMolecule, data: &str) -> Self {
        let mut ladders = HashMap::new();
        let res: Value = serde_json::from_str(data).expect("Invalid ladder JSON");
        let map = res.as_object().expect("Ladder JSON is not an object");
        for (name, parts) in map.iter() {
            ladders.insert(name.to_owned(), Ladder::new(name, parts));
        }
        Self { molecule, ladders }
    }

    pub fn molecule(&self) -> LadderMolecule {
        self.molecule
    }

    pub fn get(&self, name: &str) -> Option<&Ladder> {
        self.ladders.get(name)
    }

    pub fn iter(&self) -> impl Iterator<Item = &Ladder> {
        self.ladders.values()
    }

    pub fn names_sorted(&self) -> Vec<String> {
        let mut names = self.ladders.keys().cloned().collect::<Vec<_>>();
        names.sort_unstable();
        names
    }

    pub fn choose_for_range(
        &self,
        min_bp: usize,
        max_bp: usize,
        max_ladders: usize,
    ) -> Vec<String> {
        if self.ladders.is_empty() || max_ladders == 0 || min_bp == 0 || max_bp == 0 {
            return vec![];
        }
        let (min_bp, max_bp) = if min_bp <= max_bp {
            (min_bp, max_bp)
        } else {
            (max_bp, min_bp)
        };

        let mut ladders = self
            .iter()
            .filter_map(|ladder| {
                let lo = ladder.min_len()?;
                let hi = ladder.max_len()?;
                Some((ladder.name().to_string(), lo, hi))
            })
            .collect::<Vec<_>>();
        ladders.sort_by(|a, b| a.0.cmp(&b.0));

        let mut best_primary: Option<(String, usize, usize, usize)> = None;
        for (name, lo, hi) in &ladders {
            let uncovered = Self::range_uncovered(*lo, *hi, min_bp, max_bp);
            let span_delta = lo.abs_diff(min_bp) + hi.abs_diff(max_bp);
            let score = uncovered.saturating_mul(1000).saturating_add(span_delta);
            match &best_primary {
                None => best_primary = Some((name.clone(), *lo, *hi, score)),
                Some((_n, _l, _h, best_score)) if score < *best_score => {
                    best_primary = Some((name.clone(), *lo, *hi, score))
                }
                _ => {}
            }
        }

        let Some((primary_name, primary_lo, primary_hi, primary_score)) = best_primary else {
            return vec![];
        };
        if max_ladders == 1 {
            return vec![primary_name];
        }
        let primary_uncovered = Self::range_uncovered(primary_lo, primary_hi, min_bp, max_bp);
        if primary_uncovered == 0 {
            return vec![primary_name];
        }

        let mut best_secondary: Option<(String, usize)> = None;
        for (name, lo, hi) in ladders.iter().filter(|(n, _, _)| *n != primary_name) {
            let combined_lo = primary_lo.min(*lo);
            let combined_hi = primary_hi.max(*hi);
            let uncovered = Self::range_uncovered(combined_lo, combined_hi, min_bp, max_bp);
            let span_delta = combined_lo.abs_diff(min_bp) + combined_hi.abs_diff(max_bp);
            let score = uncovered
                .saturating_mul(1000)
                .saturating_add(span_delta)
                .saturating_add(primary_score / 10);
            match &best_secondary {
                None => best_secondary = Some((name.clone(), score)),
                Some((_n, best_score)) if score < *best_score => {
                    best_secondary = Some((name.clone(), score))
                }
                _ => {}
            }
        }

        if let Some((secondary_name, _)) = best_secondary {
            vec![primary_name, secondary_name]
        } else {
            vec![primary_name]
        }
    }

    fn range_uncovered(
        ladder_min: usize,
        ladder_max: usize,
        target_min: usize,
        target_max: usize,
    ) -> usize {
        let low_gap = target_min.saturating_sub(ladder_min);
        let high_gap = target_max.saturating_sub(ladder_max);
        low_gap.saturating_add(high_gap)
    }
}

impl Default for LadderCatalog {
    fn default() -> Self {
        default_dna_ladders()
    }
}

pub fn default_dna_ladders() -> LadderCatalog {
    LadderCatalog::from_json_str(
        LadderMolecule::Dna,
        include_str!("../assets/dna_ladders.json"),
    )
}

pub fn default_rna_ladders() -> LadderCatalog {
    LadderCatalog::from_json_str(
        LadderMolecule::Rna,
        include_str!("../assets/rna_ladders.json"),
    )
}

pub type DNALadderBand = LadderBand;
pub type DNALadder = Ladder;
pub type DNALadders = LadderCatalog;
pub type RNALadderBand = LadderBand;
pub type RNALadder = Ladder;
pub type RNALadders = LadderCatalog;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dna_ladder() {
        let dna_ladders = default_dna_ladders();
        let ladder = dna_ladders.get("GeneRuler Mix").unwrap();
        assert_eq!(ladder.bands()[1].length_bp(), 8000.0);
        assert_eq!(ladder.bands()[1].relative_strength, Some(13.0));
        assert_eq!(ladder.loading_hint(), Some(0.5));
        assert_eq!(dna_ladders.molecule(), LadderMolecule::Dna);
    }

    #[test]
    fn test_choose_for_range() {
        let dna_ladders = default_dna_ladders();
        let picks = dna_ladders.choose_for_range(120, 1300, 2);
        assert!(!picks.is_empty());
        assert!(picks.len() <= 2);
    }

    #[test]
    fn test_rna_ladder_catalog() {
        let rna_ladders = default_rna_ladders();
        assert_eq!(rna_ladders.molecule(), LadderMolecule::Rna);
        assert!(rna_ladders.get("NEB ssRNA Ladder").is_some());
        let picks = rna_ladders.choose_for_range(90, 2500, 2);
        assert!(!picks.is_empty());
        assert!(picks.len() <= 2);
    }
}
