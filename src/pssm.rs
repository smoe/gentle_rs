//use bio::stats::LogProb;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};
use std::num::ParseFloatError;
use std::f64;
//use std::path::Path;
//use bio::alphabets::dna;
//use std::fmt;

pub struct PSSM {
    pub accession: String,
    pub description: String,
    pub matrix: Vec<[f64; 4]>,
}

impl PSSM {

    pub fn new(accession: String, description: String, matrix: Vec<[f64; 4]>) -> Self {
        PSSM {
            accession,
            description,
            matrix,
        }
    }

    // Method to print the PSSM in a readable format
    pub fn print_pssm(&self) {
        println!("PSSM Accession: {}", self.accession);
        println!("PSSM Description: {}", self.description);
        let labels = ["A", "C", "G", "T"];
        for (i, row) in self.matrix.iter().enumerate() {
            println!("{}  [{:>6.2} {:>6.2} {:>6.2} {:>6.2}]", labels[i], row[0], row[1], row[2], row[3]);
        }
    }

    // Method to normalize the PSSM by its colsums
    pub fn normalize(&mut self) {
        let mut column_sums = [0.0; 4];

        // Calculate the sum of each column
        for row in &self.matrix {
            for (i, value) in row.iter().enumerate() {
                column_sums[i] += value;
            }
        }

        // Normalize each element by its column sum
        for row in &mut self.matrix {
            for (i, value) in row.iter_mut().enumerate() {
                if column_sums[i] != 0.0 {
                    *value /= column_sums[i];
                }
            }
        }
    }

    // Put Log Odds against background of 0.25
    pub fn prepare_log_odds_to_background(&mut self) {
        for row in &mut self.matrix {
            for value in row.iter_mut() {
                if *value != 0.0 {
                    *value = f64::log2(*value) + 2.0; // log2(1/4) == -2
                } else {
                    *value = -99999.99999; // substitute for -infinity
                }
            }
        }
    }
}

// Function to parse the JASPAR file and return a HashMap of PSSMs
pub fn parse_jaspar(filename: &str) -> Result<HashMap<String, PSSM>, Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let reader = io::BufReader::new(file);
    let mut pssms = HashMap::new();
    let mut lines = reader.lines();

    while let Some(line) = lines.next() {
        let line = line?;
        if line.starts_with('>') {

            // Split the header line at the first whitespace
            let mut parts = line[1..].trim().splitn(2, char::is_whitespace);
            let accession = parts.next().unwrap().to_string();
            let description = parts.next().unwrap_or("").to_string();

            // Parse the matrix
            let mut matrix = Vec::new();
            for _ in 0..4 {
                let line = lines.next().ok_or("Incomplete matrix data")??;
                let row: Result<Vec<f64>, ParseFloatError> = line
                    .split_whitespace()
                    .filter(|v| v.chars().all(|c| c.is_digit(10) || c == '.')) // Retain only numeric values
                    .map(|v| v.parse::<f64>())
                    .collect();

                match row {
                    Ok(parsed_row) if parsed_row.len() == 4 => {
                        matrix.push([parsed_row[0], parsed_row[1], parsed_row[2], parsed_row[3]]);
                    }
                    Ok(_) => {
                        return Err(Box::from("Invalid PSSM format: expected 4 numeric columns"));
                    }
                    Err(e) => {
                        return Err(Box::new(e));
                    }
                }
            }

            // Create and add PSSM to the HashMap, keyed by accession
            let pssm = PSSM::new(accession.clone(), description, matrix);
            pssms.insert(accession, pssm);

        }
    }

    Ok(pssms)
}

pub fn parse_pssm(filename: &str) -> Result<PSSM, Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let reader = io::BufReader::new(file);
    let mut lines = reader.lines();

    // Read the first line as the header and extract the name
    let name_line = lines.next().ok_or("Empty file, expected header")??;
    if !name_line.starts_with('>') {
        return Err(Box::from("Invalid PSSM format: header line does not start with '>'"));
    }

    // Split the header line at the first whitespace
    let mut parts = name_line[1..].trim().splitn(2, char::is_whitespace);
    let accession = parts.next().unwrap().to_string();
    let description = parts.next().unwrap_or("").to_string();

    // Initialize an empty matrix
    let mut matrix = Vec::new();

    // Parse the remaining lines as matrix data
    for (line_num, line) in lines.enumerate() {
        let line = line?;
        let trimmed_line = line.trim();

        // Parse each entry as a floating-point number, ignoring non-numeric characters
        let scores: Result<Vec<f64>, ParseFloatError> = trimmed_line
            .split_whitespace()
            .filter(|v| v.chars().all(|c| c.is_digit(10) || c == '.')) // Retain only numeric values
            .map(|v| v.parse::<f64>()) // Keep the Result<f64, ParseFloatError>
            .collect();

        match scores {
            Ok(parsed_scores) if parsed_scores.len() == 4 => {
                matrix.push([parsed_scores[0], parsed_scores[1], parsed_scores[2], parsed_scores[3]]);
            }
            Ok(_) => {
                eprintln!("Error: Line {} does not contain exactly 4 numeric columns: '{}'", line_num + 2, trimmed_line);
                return Err(Box::from("Invalid PSSM format: expected 4 numeric columns"));
            }
            Err(e) => {
                eprintln!("Error parsing float on line {}: '{}': {}", line_num + 2, trimmed_line, e);
                return Err(Box::new(e));
            }
        }
    }

    Ok(PSSM { accession, description, matrix })
}

pub fn score_sequence(seq: &str, pssm: &[[f64; 4]]) -> f64 {
    //let alphabet = dna::alphabet();
    let mut score = 0.0;

    for (i, base) in seq.chars().enumerate() {
        let idx = match base {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            _ => continue, // Ignore unknown bases
        };
        score += pssm[i][idx];
    }
    score
}
pub fn scan_sequence(dna_sequence: &str, pssm: &[[f64; 4]], threshold: f64) -> Vec<(usize, f64)> {
    let window_size = pssm.len();
    let mut hits = Vec::new();

    for i in 0..=dna_sequence.len() - window_size {
        let subseq = &dna_sequence[i..i + window_size];
        let score = score_sequence(subseq, pssm);

        if score >= threshold {
            hits.push((i, score));
        }
    }
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pssm_sites() -> Result<(), Box<dyn std::error::Error>> {
        let pssm_file = "test_files/MA1234.1.jaspar"; // Replace with actual file path
        let dna_sequence = "ACTGACGTACTGACGTAGCTAGCTGACGTACGTTCGATTCGA"; // Replace with actual DNA sequence
        let threshold = 0.0; // Set your desired threshold for binding site score

        // Handle the result of parse_pssm by unwrapping it
        let mut pssm = parse_pssm(pssm_file)?; // Use `?` to propagate errors if any
        println!("Parsed PSSM:");
        pssm.print_pssm(); // Print the loaded PSSM

        pssm.normalize(); // Divide values by column sums
        println!("Normalized PSSM:");
        pssm.print_pssm(); // Print the loaded PSSM

        pssm.prepare_log_odds_to_background(); // Divide values by column sums
        pssm.print_pssm(); // Print the loaded PSSM

        let hits = scan_sequence(dna_sequence, &pssm.matrix, threshold);

        for (position, score) in hits {
            println!("Potential binding site at position {}: Score = {}", position, score);
        }

        Ok(())
    }

    #[test]
    fn test_nonexistent_file() {
        let result = parse_jaspar("nonexistent_file.jaspar");
        assert!(result.is_err(), "Expected an error for a nonexistent file");

        if let Err(e) = result {
            println!("Error returned as expected: {}", e);
        }
    }

    #[test]
    fn test_jaspar_parsing() -> Result<(), Box<dyn std::error::Error>> {
        let pssm_file = "test_files/MA1234.1.jaspar"; // Replace with actual file path
        //let dna_sequence = "ACTGACGTACTGACGTAGCTAGCTGACGTACGTTCGATTCGA"; // Replace with actual DNA sequence
        //let threshold = 0.0; // Set your desired threshold for binding site score

        //let pssms = parse_jaspar("test_files/jaspar_file.jaspar")?;
        let pssms = parse_jaspar(pssm_file)?;

        for (accession, pssm) in &pssms {
            println!("\nPSSM Accession: {}", accession);
            pssm.print_pssm();
        }

        // Access a specific PSSM by name
        if let Some(pssm) = pssms.get("MA1234.1") {
            pssm.print_pssm();
        } else {
            println!("PSSM with name 'MA1234.1' not found.");
        }

        Ok(())
    }
}
