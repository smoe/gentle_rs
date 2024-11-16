use serde::Deserialize;
use std::collections::HashMap;


#[derive(Deserialize, Debug)]
pub struct PrimerConstraint {
    min_length: usize,
    max_length: usize,
    location: Option<usize>, // Optional specific start location
}

impl PrimerConstraint {
    pub fn new(min_length: usize, max_length: usize, location: Option<usize>) -> Self {
        PrimerConstraint {
            min_length,
            max_length,
            location,
        }
    }
}

/*
impl FromStr for PrimerConstraint {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut min_length = None;
        let mut max_length = None;
        let mut location = None;

        for part in s.split(',') {
            let trimmed = part.trim();
            if trimmed.starts_with("min_length=") {
                min_length = trimmed["min_length=".len()..].parse().ok();
            } else if trimmed.starts_with("max_length=") {
                max_length = trimmed["max_length=".len()..].parse().ok();
            } else if trimmed.starts_with("location=") {
                location = trimmed["location=".len()..].parse().ok();
            }
        }

        Ok(PrimerConstraint {
            min_length: min_length.ok_or("Missing min_length")?,
            max_length: max_length.ok_or("Missing max_length")?,
            location,
        })
    }
}
*/

#[derive(Deserialize, Debug)]
pub struct PrimerPairConstraint {
    forward: PrimerConstraint,
    reverse: PrimerConstraint,
    from: usize,                  // Start of the area to cover
    to: usize,                    // End of the area to cover
    max_product_size: usize,      // Maximum allowed product size
    min_temp: Option<usize>,      // Minimum annealing temperature (optional)
    max_temp: Option<usize>,      // Maximum annealing temperature (optional)
}

/*
impl FromStr for PrimerPairConstraint {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut forward_constraint = None;
        let mut reverse_constraint = None;
        let mut from = None;
        let mut to = None;
        let mut max_product_size = None;
        let mut min_temp = None;
        let mut max_temp = None;

        for part in s.split(',') {
            let trimmed = part.trim();

            if trimmed.starts_with("forward") {
                let fwd_text = trimmed["forward=".len()..].trim();
                forward_constraint = Some(fwd_text.parse::<PrimerConstraint>().map_err(|e| e.to_string())?);
            } else if trimmed.starts_with("reverse") {
                let rev_text = trimmed["reverse=".len()..].trim();
                reverse_constraint = Some(rev_text.parse::<PrimerConstraint>().map_err(|e| e.to_string())?);
            } else if trimmed.starts_with("from=") {
                from = trimmed["from=".len()..].parse().ok();
            } else if trimmed.starts_with("to=") {
                to = trimmed["to=".len()..].parse().ok();
            } else if trimmed.starts_with("max_product_size=") {
                max_product_size = trimmed["max_product_size=".len()..].parse().ok();
            } else if trimmed.starts_with("min_temp=") {
                min_temp = trimmed["min_temp=".len()..].parse().ok();
            } else if trimmed.starts_with("max_temp=") {
                max_temp = trimmed["max_temp=".len()..].parse().ok();
            }
        }

        Ok(PrimerPairConstraint {
            forward_constraint: forward_constraint.ok_or("Missing forward primer constraint")?,
            reverse_constraint: reverse_constraint.ok_or("Missing reverse primer constraint")?,
            from: from.ok_or("Missing from constraint")?,
            to: to.ok_or("Missing to constraint")?,
            max_product_size: max_product_size.ok_or("Missing max_product_size constraint")?,
            min_temp,
            max_temp,
        })
    }
}
*/

pub fn count_occurrences(sequence: &str, primer: &str) -> usize {
    sequence.matches(primer).count()
}

// Function to estimate "melting temperature" based on GC content
pub fn estimate_melting_temp(primer: &str) -> usize {
    primer.chars().filter(|&c| c == 'G' || c == 'C').count()
}

// Reverse complement function
pub fn reverse_complement(dna: &str) -> String {
    dna.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => c,
        })
        .collect()
}

// Generate primers for a specific constraint
pub fn generate_primers(sequence: &str, constraint: &PrimerConstraint) -> Vec<String> {
    let mut primers = Vec::new();
    for start in 0..sequence.len() {
        if let Some(location) = constraint.location {
            if start != location {
                continue;
            }
        }
        for len in constraint.min_length..=constraint.max_length {
            if start + len <= sequence.len() {
                primers.push(sequence[start..start + len].to_string());
            }
        }
    }
    primers
}

// Generate valid primer pairs based on pair constraints
pub fn generate_primer_pairs(sequence: &str, pair_constraint: &PrimerPairConstraint) -> Vec<(String, String)> {
    let mut primer_pairs = Vec::new();

    // Generate the reverse complement of the entire sequence
    let reverse_sequence = reverse_complement(sequence);

    let forward_primers = generate_primers(sequence, &pair_constraint.forward);
    let reverse_primers = generate_primers(sequence, &pair_constraint.reverse)
        .into_iter()
        .map(|s| reverse_complement(&s))  // Get reverse complement of each reverse primer
        .collect::<Vec<_>>();

    for forward_primer in &forward_primers {
        let forward_len = forward_primer.len();
        let forward_temp = estimate_melting_temp(forward_primer);

        // Ensure forward primer is unique in the sequence
        if count_occurrences(sequence, forward_primer) != 1 {
            continue;
        }

        for reverse_primer in &reverse_primers {
            let reverse_len = reverse_primer.len();
            let reverse_temp = estimate_melting_temp(reverse_primer);

            // Ensure reverse primer is unique in the reverse complement sequence
            if count_occurrences(&reverse_sequence, reverse_primer) != 1 {
                continue;
            }

            // Check primer length difference constraint (within 2 nucleotides)
            if (forward_len as isize - reverse_len as isize).abs() > 2 {
                continue;
            }

            // Check approximate melting temperature similarity (within 1 GC count)
            if (forward_temp as isize - reverse_temp as isize).abs() > 1 {
                continue;
            }

            // Check optional min and max temperature constraints
            if let Some(min_temp) = pair_constraint.min_temp {
                if forward_temp < min_temp || reverse_temp < min_temp {
                    continue;
                }
            }
            if let Some(max_temp) = pair_constraint.max_temp {
                if forward_temp > max_temp || reverse_temp > max_temp {
                    continue;
                }
            }

            // Ensure forward and reverse primers meet product size constraint
            let fwd_start = sequence.find(forward_primer).unwrap();
            let rev_start = reverse_sequence.find(reverse_primer).unwrap();
            let product_size = (rev_start + reverse_len) + fwd_start;

            if product_size <= pair_constraint.max_product_size {
                primer_pairs.push((forward_primer.clone(), reverse_primer.clone()));
            }
        }
    }

    primer_pairs
}

/*
 * Cleavage
 */


#[derive(Deserialize, Debug)]
struct CleavageConstraint {
    enzyme: String,           // Enzyme name (e.g., "EcoRI")
    site: String,             // Recognition site (e.g., "GAATTC")
    allow: bool               // Whether to allow cleavage at this site
}

#[derive(Deserialize, Debug)]
struct CleavageSiteConstraint {
    enzymes: Vec<CleavageConstraint>,  // List of enzyme-specific constraints
    from: Option<usize>,               // Start of region to consider
    to: Option<usize>,                 // End of region to consider
}

// Function to find cleavage sites in the sequence based on the given constraints
fn find_cleavage_sites(sequence: &str, constraint: &CleavageSiteConstraint) -> HashMap<String, Vec<usize>> {
    let mut cleavage_sites: HashMap<String, Vec<usize>> = HashMap::new();

    for enzyme in &constraint.enzymes {
        let mut sites = Vec::new();
        let mut start = 0;

        while let Some(pos) = sequence[start..].find(&enzyme.site) {
            let abs_pos = start + pos;
            start = abs_pos + 1;

            // Only add the site if allowed
            if enzyme.allow {
                sites.push(abs_pos);
            }
        }

        if !sites.is_empty() {
            cleavage_sites.insert(enzyme.enzyme.clone(), sites);
        }
    }

    cleavage_sites
}

// Function to find the closest enzyme sites to a region of interest (ROI)
fn closest_cleavage_sites(sequence: &str, constraint: &CleavageSiteConstraint, roi_start: usize, roi_end: usize) -> HashMap<String, Option<(usize, usize)>> {
    let cleavage_sites = find_cleavage_sites(sequence, constraint);
    let mut closest_sites: HashMap<String, Option<(usize, usize)>> = HashMap::new();

    for (enzyme, sites) in cleavage_sites {
        let mut closest_start = None;
        let mut closest_end = None;
        let mut min_start_dist = usize::MAX;
        let mut min_end_dist = usize::MAX;

        for &site in &sites {
            let start_dist = if site <= roi_start { roi_start - site } else { site - roi_start };
            let end_dist = if site <= roi_end { roi_end - site } else { site - roi_end };

            if start_dist < min_start_dist {
                min_start_dist = start_dist;
                closest_start = Some(site);
            }
            if end_dist < min_end_dist {
                min_end_dist = end_dist;
                closest_end = Some(site);
            }
        }

        closest_sites.insert(enzyme, Some((closest_start.unwrap_or(roi_start), closest_end.unwrap_or(roi_end))));
    }

    closest_sites
}




/*
 * Tests
 */


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constraints() -> Result<(), Box<dyn std::error::Error>> {
        // Sample DNA sequence
        let sequence = "ATGCGTACGTTAGCTAGCTTACGGTAGCTAG";
        let primer_pair_constraint_json = r#"
{
    "forward": {
        "min_length": 5,
        "max_length": 10,
        "location": 3
    },
    "reverse": {
        "min_length": 5,
        "max_length": 10
    },
    "from": 3,
    "to": 20,
    "max_product_size": 100,
    "min_temp": 3,
    "max_temp": 7
}
"#;        
        // Parse JSON into a PrimerPairConstraint struct
        let pair_constraint: PrimerPairConstraint = serde_json::from_str(primer_pair_constraint_json)?;
        // Print parsed constraints for verification
        println!("{:?}", pair_constraint);

        // Generate primer pairs based on constraints
        let primer_pairs = generate_primer_pairs(sequence, &pair_constraint);

        println!("Generated primer pairs:");
        for (forward, reverse) in &primer_pairs {
            println!("Forward: {}, Reverse: {}", forward, reverse);
        }

        // Test primer pair generation
        assert!(!primer_pairs.is_empty());
        Ok(())
    }

    fn test_helper_functions() {
        let sequence = "ATGCGTACGTTAGC";
        // Test individual functions
        assert_eq!(count_occurrences(sequence, "ATGC"), 1);
        assert_eq!(estimate_melting_temp("GCGT"), 3);
        assert_eq!(reverse_complement("ATGC"), "GCAT");

        // Test primer generation
        let primers = generate_primers(sequence, &PrimerConstraint::new(5, 10, None));
        assert!(!primers.is_empty());

    }

    // Cleavage

    #[test]
    fn test_cleavage_constraints() -> Result<(), Box<dyn std::error::Error>> {
        // Sample DNA sequence
        let sequence = "ATGCGTACGTTAGCTAGAATTCTAGTACCGGATCCAGCTAGTGAATTCGTAGGATCC";

        // JSON input defining enzymes and constraints
        let cleavage_constraint_json = r#"
        {
            "enzymes": [
                {
                    "enzyme": "EcoRI",
                    "site": "GAATTC",
                    "allow": true
                },
                {
                    "enzyme": "BamHI",
                    "site": "GGATCC",
                    "allow": true
                }
            ],
            "from": 0,
            "to": 60
        }
        "#;

        // Parse JSON into CleavageSiteConstraint struct
        let constraint: CleavageSiteConstraint = serde_json::from_str(cleavage_constraint_json)?;

        // Define the Region of Interest (ROI)
        let roi_start = 20;
        let roi_end = 40;

        // Find the closest cleavage sites to the ROI
        let closest_sites = closest_cleavage_sites(sequence, &constraint, roi_start, roi_end);

        // Print results for verification
        println!("Closest cleavage sites to the ROI ({}, {}):", roi_start, roi_end);
        for (enzyme, sites) in closest_sites {
            if let Some((closest_start, closest_end)) = sites {
                println!("Enzyme: {}, Closest Start: {}, Closest End: {}", enzyme, closest_start, closest_end);
            } else {
                println!("Enzyme: {}, No cleavage sites found near ROI", enzyme);
            }
        }

        Ok(())
    }
}


