use anyhow::{anyhow, Result};
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, Seek, SeekFrom, Write};
use std::num::ParseFloatError;
use serde_json::Value;

// Should consider falling back on https://docs.rs/bio/latest/src/bio/pattern_matching/pssm/mod.rs.html

#[derive(Debug, Clone, Default)]
pub struct PSSM {
    accession: String,
    description: String,
    matrix: Vec<Vec<f64>>,
    is_normalized: bool,
    log_odds_prepared: bool,
}

impl PSSM {
    pub fn new(accession: String, description: String, matrix: Vec<Vec<f64>>) -> Self {
        PSSM {
            accession,
            description,
            matrix,
            ..Default::default()
        }
    }

    pub fn accession(&self) -> &str {
        &self.accession
    }

    pub fn description(&self) -> &str {
        &self.description
    }

    // Method to normalize the PSSM by its colsums
    pub fn normalize(&mut self) {
        if self.is_normalized {
            return;
        }
        self.is_normalized = true;
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
        if self.log_odds_prepared {
            return;
        }
        self.log_odds_prepared = true;
        for row in &mut self.matrix {
            for value in row.iter_mut() {
                if (*value) != 0.0 {
                    *value = f64::log2(*value) + 2.0; // log2(1/4) == -2
                } else {
                    *value = f64::NEG_INFINITY;
                }
            }
        }
    }

    pub fn from_pssm_file(filename: &str) -> Result<Self> {
        let file = File::open(filename)?;
        let reader = io::BufReader::new(file);
        let mut lines = reader.lines();

        // Read the first line as the header and extract the name
        let name_line = lines
            .next()
            .ok_or(anyhow!("Empty file, expected header"))??;
        if !name_line.starts_with('>') {
            return Err(anyhow!(
                "Invalid PSSM format: header line does not start with '>'",
            ));
        }

        // Split the header line at the first whitespace
        let mut parts = name_line[1..].trim().splitn(2, char::is_whitespace);
        let accession = parts.next().unwrap().to_string();
        let description = parts.next().unwrap_or("").to_string();

        // Initialize an empty matrix
        let mut matrix = Vec::new();

        // Parse the remaining lines as matrix data
        for (line_num, line) in lines.enumerate() {
            Self::parse_pssm_line(line, &mut matrix, line_num)?;
        }

        Ok(PSSM {
            accession,
            description,
            matrix,
            ..Default::default()
        })
    }

    pub fn from_elixir_api(id: &str) -> Result<HashMap<String, PSSM>> {
        let url = format!("https://jaspar.elixir.no/api/v1/matrix/{id}/?format=jaspar");
        let text = reqwest::blocking::get(url)?.text()?;
        let mut file = tempfile::tempfile()?;
        file.write_all(text.as_bytes())?;
        file.seek(SeekFrom::Start(0))?;

        let reader = io::BufReader::new(file);
        let mut pssms = HashMap::new();
        let mut lines = reader.lines();

        while let Some(line) = lines.next() {
            Self::parse_jaspar_line(line, &mut lines, &mut pssms)?;
        }

        Ok(pssms)
    }

    // Function to parse the JASPAR file and return a HashMap of PSSMs
    pub fn from_jaspar_file(filename: &str) -> Result<HashMap<String, PSSM>> {
        let file = File::open(filename)?;
        let reader = io::BufReader::new(file);
        let mut pssms = HashMap::new();
        let mut lines = reader.lines();

        while let Some(line) = lines.next() {
            Self::parse_jaspar_line(line, &mut lines, &mut pssms)?;
        }

        Ok(pssms)
    }

    pub fn from_json_file(filename: &str) -> Result<HashMap<String, PSSM>> {
        let file = File::open(filename)?;
        let reader = io::BufReader::new(file);
        let json: Value = serde_json::from_reader(reader)?;

        let mut pssms = HashMap::new();

        if let Some(entries) = json.as_array() {
            for entry in entries {
                let accession = entry["matrix_id"].as_str().ok_or(anyhow!("Missing matrix_id"))?.to_string();
                let description = entry["name"].as_str().unwrap_or("").to_string();
                let mut matrix = vec![vec![]; 4];

                for (base, counts) in entry["pfm"].as_object().ok_or(anyhow!("Invalid pfm format"))? {
                    let row = match base.as_str() {
                        "A" => &mut matrix[0],
                        "C" => &mut matrix[1],
                        "G" => &mut matrix[2],
                        "T" => &mut matrix[3],
                        _ => return Err(anyhow!("Invalid base: {}", base)),
                    };

                    for count in counts.as_array().ok_or(anyhow!("Invalid counts format"))? {
                        let count_value = count.as_f64().or_else(|| count.as_str().and_then(|s| s.parse::<f64>().ok()))
                            .ok_or(anyhow!("Invalid count value: {:?}", count))?;
                        row.push(count_value);
                    }
                }

                let pssm = PSSM::new(accession.clone(), description, matrix);
                pssms.insert(accession, pssm);
            }
        }

        Ok(pssms)
    }

    pub fn scan_sequence(&self, dna_sequence: &str, threshold: f64) -> Vec<(usize, f64)> {
        let window_size = self.matrix.len();
        let mut hits = Vec::new();

        for i in 0..=dna_sequence.len() - window_size {
            let subseq = &dna_sequence[i..i + window_size];
            let score = self.score_sequence(subseq);

            if score >= threshold {
                hits.push((i, score));
            }
        }
        hits
    }

    fn score_sequence(&self, seq: &str) -> f64 {
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
            score += self.matrix[i][idx];
        }
        score
    }

    fn parse_pssm_line(
        line: std::result::Result<String, io::Error>,
        matrix: &mut Vec<Vec<f64>>,
        line_num: usize,
    ) -> Result<()> {
        let line = line?;
        let trimmed_line = line.trim();
        let scores: Result<Vec<f64>, ParseFloatError> = trimmed_line
            .split_whitespace()
            .filter(|v| v.chars().all(|c| c.is_ascii_digit() || c == '.')) // Retain only numeric values
            .map(|v| v.parse::<f64>()) // Keep the Result<f64, ParseFloatError>
            .collect();
        match scores {
            Ok(parsed_scores) if parsed_scores.len() == 4 => {
                matrix.push(parsed_scores);
            }
            Ok(_) => {
                eprintln!(
                    "Error: Line {} does not contain exactly 4 numeric columns: '{}'",
                    line_num + 2,
                    trimmed_line
                );
                return Err(anyhow!("Invalid PSSM format: expected 4 numeric columns"));
            }
            Err(e) => {
                eprintln!(
                    "Error parsing float on line {}: '{}': {}",
                    line_num + 2,
                    trimmed_line,
                    e
                );
                return Err(anyhow!(e));
            }
        }
        Ok(())
    }

    fn parse_jaspar_line(
        line: std::result::Result<String, io::Error>,
        lines: &mut io::Lines<io::BufReader<File>>,
        pssms: &mut HashMap<String, PSSM>,
    ) -> Result<()> {
        let line = line?;
        if let Some(end) = line.strip_prefix('>') {
            // Split the header line at the first whitespace
            let mut parts = end.trim().splitn(2, char::is_whitespace);
            let accession = parts.next().unwrap().to_string();
            let description = parts.next().unwrap_or("").to_string();

            // Parse the matrix
            let mut matrix = Vec::new();
            for _ in 0..4 {
                let line = lines.next().ok_or(anyhow!("Incomplete matrix data"))??;
                let row: Result<Vec<f64>, ParseFloatError> = line
                    .split_whitespace()
                    .filter(|v| v.chars().all(|c| c.is_ascii_digit() || c == '.')) // Retain only numeric values
                    .map(|v| v.parse::<f64>())
                    .collect();

                match row {
                    Ok(parsed_row) => {
                        matrix.push(parsed_row);
                    }
                    Err(e) => {
                        return Err(anyhow!(e));
                    }
                }
            }

            // Create and add PSSM to the HashMap, keyed by accession
            let pssm = PSSM::new(accession.clone(), description, matrix);
            pssms.insert(accession, pssm);
        }
        Ok(())
    }
}

impl fmt::Display for PSSM {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "PSSM Accession: {}", self.accession)?;
        writeln!(f, "PSSM Description: {}", self.description)?;
        let labels = ["A", "C", "G", "T"];
        for (i, row) in self.matrix.iter().enumerate() {
            let display_row = row
                .iter()
                .map(|v| format!("{:>6.2}", v))
                .collect::<Vec<String>>()
                .join(" ");
            writeln!(f, "{}  [{display_row}]", labels[i])?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pssm_sites() -> Result<()> {
        let pssm_file = "test_files/MA1234.1.jaspar"; // JASPAR file
        let dna_sequence = "ACTGACGTACTGACGTAGCTAGCTGACGTACGTTCGATTCGA"; // Replace with actual DNA sequence
        let threshold = 0.0; // Set your desired threshold for binding site score

        // Handle the result of PSSM::from_file by unwrapping it
        let mut pssm = PSSM::from_pssm_file(pssm_file)?; // Use `?` to propagate errors if any
        assert_eq!(pssm.accession, "MA1234.1");
        assert_eq!(pssm.description, "ACGT");
        assert_eq!(pssm.matrix[1][1], 4429.0);

        // Divide values by column sums
        pssm.normalize();
        assert_eq!(pssm.matrix[1][1], 1.0);

        // Divide values by column sums
        pssm.prepare_log_odds_to_background();
        assert_eq!(pssm.matrix[1][1], 2.0);
        assert_eq!(pssm.matrix[0][1], f64::NEG_INFINITY);

        let hits = pssm.scan_sequence(dna_sequence, threshold);
        assert_eq!(hits, [(4, 8.0), (12, 8.0), (25, 8.0), (29, 8.0)]);

        Ok(())
    }

    #[test]
    fn test_nonexistent_file() {
        let result = PSSM::from_jaspar_file("nonexistent_file.jaspar");
        assert!(result.is_err(), "Expected an error for a nonexistent file");
    }

    #[test]
    fn test_elixir_api() {
        let pssms = PSSM::from_elixir_api("MA0265.1").unwrap();
        assert_eq!(pssms.len(), 1);
        let pssm = pssms.get("MA0265.1").unwrap();
        assert_eq!(pssm.accession, "MA0265.1");
        assert_eq!(pssm.description, "ABF1");
        assert_eq!(pssm.matrix[1][5], 24.0);
    }

    #[test]
    fn test_jaspar_parsing() -> Result<()> {
        let pssm_file = "test_files/MA1234.1.jaspar"; // Replace with actual file path
        let pssms = PSSM::from_jaspar_file(pssm_file)?;
        assert_eq!(pssms.len(), 1);

        // Access a specific PSSM by name
        let pssm = pssms.get("MA1234.1").unwrap();

        // Check indivudual value
        assert_eq!(pssm.matrix[1][1], 4429.00);

        Ok(())
    }

    #[test]
    fn test_json_parsing() -> Result<()> {
        let json_file = "assets/jaspar_2022.json"; // Replace with actual file path
        let pssms = PSSM::from_json_file(json_file)?;
        
        assert_eq!(pssms.len(), 1956); // Equals `grep -c "^>" data/JASPAR_2022.txt`

        // Access a specific PSSM by name
        let pssm = pssms.get("MA0004.1").unwrap();

        // Check individual value
        assert_eq!(pssm.matrix[1][1], 0.0);
        assert_eq!(pssm.matrix[2][3], 20.0);

        let pssm = pssms.get("MA0006.1").unwrap();
        assert_eq!(pssm.matrix[0][0], 3.0);
        assert_eq!(pssm.matrix[1][1], 0.0);
        assert_eq!(pssm.matrix[0][2], 0.0);
        assert_eq!(pssm.matrix[3][0], 11.0);
        assert_eq!(pssm.matrix[1][2], 23.0);
        assert_eq!(pssm.matrix[2][5], 24.0);

        Ok(())
    }
}
