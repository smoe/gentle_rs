use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Write, BufReader, BufRead}; // Hinzugefügt: BufRead
use reqwest::blocking::get;
use std::path::Path;
use flate2::read::GzDecoder;

/// Management of locally available (typically public) genome data

#[derive(Default, Deserialize, Serialize, Debug, Clone)]
#[allow(dead_code)]
struct GenomeDataLocal {
    ncbi_taxonomy_id: Option<u32>,
    description: Option<String>,
    sequence_remote: Option<String>,
    annotations_remote: Option<String>,
    sequence_local: Option<String>,
    annotations_local: Option<String>,
    #[serde(default = "default_cache_dir")]
    cache_dir: Option<String>,
    #[serde(default)]
    found_locally: bool,
    #[serde(default)]
    downloaded: bool,
}

fn default_cache_dir() -> Option<String> {
    Some("/tmp".to_string())
}

#[derive(Debug, Deserialize, Serialize, Clone)]
#[serde(rename_all = "snake_case")]
enum GeneBiotype {
    #[serde(rename = "undefined")]
    Undefined,
    NcRNA,
    ProteinCoding,
    Pseudogene,
    RRNA,
    SnRNA,
    SnoRNA,
    TRNA,
    TransposableElement,
    Other,
}

impl Default for GeneBiotype {
    fn default() -> Self {
        GeneBiotype::Undefined
    }
}

#[derive(Debug, Deserialize, Serialize, Clone)]
#[serde(rename_all = "snake_case")]
enum TranscriptBiotype {
    #[serde(rename = "undefined")]
    Undefined,
    NcRNA,
    ProteinCoding,
    Pseudogene,
    RRNA,
    SnRNA,
    SnoRNA,
    TRNA,
    TransposableElement,
    Other,
}

impl Default for TranscriptBiotype {
    fn default() -> Self {
        TranscriptBiotype::Undefined
    }
}

#[derive(Debug, Default, Clone)]
#[allow(dead_code)]
struct Exon {
    start: usize,
    end: usize,
    exon_number: Option<String>,
    exon_id: Option<String>,
}

#[derive(Debug, Default, Clone)]
#[allow(dead_code)]
struct Transcript {
    transcript_id: String,
    exons: Vec<Exon>,
    protein_id: Option<String>,
    transcript_biotype: TranscriptBiotype,
}

#[derive(Debug, Default, Clone)]
#[allow(dead_code)]
struct Gene {
    gene_id: String,
    gene_name: Option<String>,
    chromosome: String,
    start: usize,
    end: usize,
    gene_biotype: GeneBiotype,
    transcripts: Vec<Transcript>,
    exons: Vec<Exon>,
}

#[derive(Debug, Default, Clone)]
#[allow(dead_code)]
struct Protein {
    protein_id: String,
    gene_id: String,
    transcript_id: String,
}

#[derive(Debug, Default)]
#[allow(dead_code)]
struct GenomeContainer {
    genomes: HashMap<String, HashMap<String, String>>, // Genome -> Chromosome -> Sequence
    genomes_metadata: HashMap<String, GenomeDataLocal>, // Genome -> GenomeDataLocal
    genes: HashMap<String, Gene>, // Gene ID -> Gene
    transcripts: HashMap<String, Transcript>, // Transcript ID -> Transcript
    proteins: HashMap<String, Protein>, // Protein ID -> Protein
}

#[allow(dead_code)]
impl GenomeContainer {
    fn new() -> Self {
        GenomeContainer {
            genomes: HashMap::new(),
            genomes_metadata: HashMap::new(),
            genes: HashMap::new(),
            transcripts: HashMap::new(),
            proteins: HashMap::new(),
        }
    }

    fn from_json_file(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let data = fs::read_to_string(path)?;
        let mut genomes_metadata: HashMap<String, GenomeDataLocal> = serde_json::from_str(&data)?;

        for metadata in genomes_metadata.values_mut() {
            let sequence_exists = metadata.sequence_local.as_ref().map_or(false, |path| Path::new(path).exists());
            let annotations_exists = metadata.annotations_local.as_ref().map_or(false, |path| Path::new(path).exists());

            if sequence_exists && annotations_exists {
                metadata.downloaded = true;
                metadata.found_locally = true;
            } else if sequence_exists || annotations_exists {
                return Err(Box::new(io::Error::new(io::ErrorKind::Other, "Only one of the files exists locally")));
            }
        }

        let mut container = GenomeContainer::new();
        for (genome_name, metadata) in genomes_metadata {
            container.add_genome(genome_name, metadata)?;
        }
        Ok(container)
    }

    fn add_genome(&mut self, genome_name: String, genomes_metadata: GenomeDataLocal) -> Result<(), String> {
        if self.genomes.contains_key(&genome_name) || self.genomes_metadata.contains_key(&genome_name) {
            return Err(format!("Genome '{}' already exists", genome_name));
        }
        self.genomes_metadata.insert(genome_name.clone(), genomes_metadata);
        self.genomes.insert(genome_name, HashMap::new());
        Ok(())
    }

    fn add_chromosome(&mut self, genome_name: &str, chromosome_name: String, sequence: String) -> Result<(), String> {
        if let Some(chromosomes) = self.genomes.get_mut(genome_name) {
            chromosomes.insert(chromosome_name, sequence);
            Ok(())
        } else {
            Err(format!("Genome '{}' not found", genome_name))
        }
    }

    fn get_sequence(&self, genome_name: &str, chromosome_name: &str) -> Result<&String, String> {
        self.genomes.get(genome_name)
            .ok_or_else(|| format!("Genome '{}' not found", genome_name))?
            .get(chromosome_name)
            .ok_or_else(|| format!("Chromosome '{}' not found in genome '{}'", chromosome_name, genome_name))
    }

    fn get_genomes_metadata(&self, genome_name: &str) -> Result<&GenomeDataLocal, String> {
        self.genomes_metadata.get(genome_name)
            .ok_or_else(|| format!("Genome '{}' not found", genome_name))
    }

    fn download_genome_files(&mut self, genome_name: &str) -> Result<(), String> {
        let metadata = self.genomes_metadata.get_mut(genome_name)
            .ok_or_else(|| format!("Genome '{}' not found", genome_name))?;

        if metadata.found_locally {
            return Ok(());
        }

        if metadata.downloaded {
            return Err("Genome already downloaded".to_string());
        }

        let mut sequence_exists = false;
        let mut annotations_exists = false;

        if let Some(sequence_url) = &metadata.sequence_remote {
            let file_name = Path::new(sequence_url)
                .file_name()
                .ok_or("Invalid sequence URL")?
                .to_str()
                .ok_or("Invalid sequence URL")?;
            let local_path = metadata.sequence_local.get_or_insert_with(|| format!("/tmp/{}", file_name));
            if Path::new(local_path).exists() {
                sequence_exists = true;
            } else {
                let response = get(sequence_url).map_err(|e| e.to_string())?;
                let mut file = File::create(local_path).map_err(|e| e.to_string())?;
                file.write_all(&response.bytes().map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
            }
        }

        if let Some(annotations_url) = &metadata.annotations_remote {
            let file_name = Path::new(annotations_url)
                .file_name()
                .ok_or("Invalid annotations URL")?
                .to_str()
                .ok_or("Invalid annotations URL")?;
            let local_path = metadata.annotations_local.get_or_insert_with(|| format!("/tmp/{}", file_name));
            if Path::new(local_path).exists() {
                annotations_exists = true;
            } else if !sequence_exists {
                let response = get(annotations_url).map_err(|e| e.to_string())?;
                let mut file = File::create(local_path).map_err(|e| e.to_string())?;
                file.write_all(&response.bytes().map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
            }
        }

        if sequence_exists && annotations_exists {
            metadata.found_locally = true;
        } else if sequence_exists || annotations_exists {
            return Err("Only one of the files exists locally".to_string());
        }

        metadata.downloaded = true;
        Ok(())
    }

    fn load_fasta_sequences(&mut self, genome_name: &str) -> Result<(), String> {
        let metadata = self.genomes_metadata.get(genome_name)
            .ok_or_else(|| format!("Genome '{}' not found", genome_name))?;

        let sequence_path = metadata.sequence_local.as_ref().ok_or("Local sequence path not set")?;
        let file = File::open(sequence_path).map_err(|e| e.to_string())?;
        let reader = BufReader::new(GzDecoder::new(file));

        let mut current_chromosome = String::new();
        let mut current_sequence = String::new();

        for line in reader.lines() {
            let line = line.map_err(|e| e.to_string())?;
            if line.starts_with('>') {
                if !current_chromosome.is_empty() {
                    self.add_chromosome(genome_name, current_chromosome.clone(), current_sequence.clone())?;
                }
                current_chromosome = line[1..].split_whitespace().next().unwrap().to_string();
                current_sequence.clear();
            } else {
                current_sequence.push_str(&line);
            }
        }

        if !current_chromosome.is_empty() {
            self.add_chromosome(genome_name, current_chromosome, current_sequence)?;
        }

        Ok(())
    }

    fn load_gff_annotations(&mut self, genome_name: &str) -> Result<(), String> {
        let metadata = self.genomes_metadata.get(genome_name)
            .ok_or_else(|| format!("Genome '{}' not found", genome_name))?;

        let annotations_path = metadata.annotations_local.as_ref().ok_or("Local annotations path not set")?;
        let file = File::open(annotations_path).map_err(|e| e.to_string())?;
        let reader = BufReader::new(GzDecoder::new(file));

        for line in reader.lines() {
            let line = line.map_err(|e| e.to_string())?;
            if line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != 9 {
                return Err("Invalid GFF format - not experiencing the 9 tab-separated columns".to_string());
            }

            let chromosome = fields[0].to_string();
            let feature_type = fields[2].to_string();
            let start: usize = fields[3].parse().map_err(|e: std::num::ParseIntError| e.to_string())?;
            let end: usize = fields[4].parse().map_err(|e: std::num::ParseIntError| e.to_string())?;
            let attributes = fields[8].to_string();

            let mut gene_id = None;
            let mut transcript_id = None;
            let mut exon_number = None;
            let mut exon_id = None;
            let mut protein_id = None;
            let mut gene_name = None;
            let mut gene_biotype = None;
            let mut transcript_biotype = None;

            for attribute in attributes.split(';') {
                let parts: Vec<&str> = attribute.trim().split_whitespace().collect();
                if parts.len() == 2 {
                    match parts[0] {
                        "gene_id" => gene_id = Some(parts[1].replace("\"", "")),
                        "transcript_id" => transcript_id = Some(parts[1].replace("\"", "")),
                        "exon_number" => exon_number = Some(parts[1].replace("\"", "")),
                        "exon_id" => exon_id = Some(parts[1].replace("\"", "")),
                        "protein_id" => protein_id = Some(parts[1].replace("\"", "")),
                        "gene_name" => gene_name = Some(parts[1].replace("\"", "")),
                        "gene_biotype" => gene_biotype = Some(parts[1].replace("\"", "")),
                        "transcript_biotype" => transcript_biotype = Some(parts[1].replace("\"", "")),
                        _ => {}
                    }
                }
            }

            let gene_biotype_enum = match gene_biotype.as_deref() {
                Some("ncRNA") => GeneBiotype::NcRNA,
                Some("protein_coding") => GeneBiotype::ProteinCoding,
                Some("pseudogene") => GeneBiotype::Pseudogene,
                Some("rRNA") => GeneBiotype::RRNA,
                Some("snRNA") => GeneBiotype::SnRNA,
                Some("snoRNA") => GeneBiotype::SnoRNA,
                Some("tRNA") => GeneBiotype::TRNA,
                Some("transposable_element") => GeneBiotype::TransposableElement,
                _ => GeneBiotype::Other,
            };

            let transcript_biotype_enum = match transcript_biotype.as_deref() {
                Some("ncRNA") => TranscriptBiotype::NcRNA,
                Some("protein_coding") => TranscriptBiotype::ProteinCoding,
                Some("pseudogene") => TranscriptBiotype::Pseudogene,
                Some("rRNA") => TranscriptBiotype::RRNA,
                Some("snRNA") => TranscriptBiotype::SnRNA,
                Some("snoRNA") => TranscriptBiotype::SnoRNA,
                Some("tRNA") => TranscriptBiotype::TRNA,
                Some("transposable_element") => TranscriptBiotype::TransposableElement,
                _ => TranscriptBiotype::Other,
            };

            match feature_type.as_str() {
                "gene" => {
                    if let Some(gene_id) = gene_id {
                        self.genes.insert(gene_id.clone(), Gene {
                            gene_id,
                            gene_name,
                            chromosome: chromosome.clone(),
                            start,
                            end,
                            gene_biotype: gene_biotype_enum,
                            transcripts: Vec::new(),
                            exons: Vec::new(),
                        });
                    }
                }
                "transcript" => {
                    if let (Some(gene_id), Some(transcript_id)) = (gene_id.clone(), transcript_id.clone()) {
                        let transcript = Transcript {
                            transcript_id: transcript_id.clone(),
                            exons: Vec::new(),
                            protein_id: protein_id.clone(),
                            transcript_biotype: transcript_biotype_enum,
                        };
                        self.transcripts.insert(transcript_id.clone(), transcript.clone());
                        if let Some(gene) = self.genes.get_mut(&gene_id) {
                            gene.transcripts.push(transcript);
                            if let Some(protein_id) = protein_id {
                                self.proteins.insert(protein_id.clone(), Protein {
                                    protein_id,
                                    gene_id: gene_id.clone(),
                                    transcript_id: transcript_id.clone(),
                                });
                            }
                        }
                    }
                }
                "exon" => {
                    if let (Some(gene_id), Some(transcript_id)) = (gene_id.clone(), transcript_id.clone()) {
                        if let Some(gene) = self.genes.get_mut(&gene_id) {
                            let exon = Exon {
                                start,
                                end,
                                exon_number: exon_number.clone(),
                                exon_id: exon_id.clone(),
                            };
                            if let Some(transcript) = self.transcripts.get_mut(&transcript_id) {
                                transcript.exons.push(exon.clone());
                            }
                            gene.exons.push(exon);
                        }
                    }
                }
                _ => {}
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_json_file() -> Result<(), Box<dyn std::error::Error>> {
        let config_path = "assets/genomes.json";
        let genomes_metadata = GenomeContainer::from_json_file(config_path)?;
        println!("Loaded config: {:?}", genomes_metadata);
        Ok(())
    }

    #[test]
    fn test_genome_container() -> Result<(), Box<dyn std::error::Error>> {
        let config_path = "assets/genomes.json";
        let mut container = GenomeContainer::from_json_file(config_path)?;

        container.add_chromosome("Human GRCh38 Ensembl 113", "1".to_string(), "ATCG".to_string())?;
        container.add_chromosome("Human GRCh38 Ensembl 113", "2".to_string(), "GCTA".to_string())?;
        container.add_chromosome("Human GRCh38 Ensembl 113", "X".to_string(), "CGTA".to_string())?;
        container.add_chromosome("Human GRCh38 Ensembl 113", "Y".to_string(), "TACG".to_string())?;

        assert_eq!(container.get_sequence("Human GRCh38 Ensembl 113", "1")?, "ATCG");
        assert_eq!(container.get_sequence("Human GRCh38 Ensembl 113", "2")?, "GCTA");
        assert_eq!(container.get_sequence("Human GRCh38 Ensembl 113", "X")?, "CGTA");
        assert_eq!(container.get_sequence("Human GRCh38 Ensembl 113", "Y")?, "TACG");
        assert!(container.get_genomes_metadata("Human GRCh38 Ensembl 113").is_ok());

        Ok(())
    }

    #[test]
    fn test_add_chromosome_to_nonexistent_genome() {
        let mut container = GenomeContainer::new();
        let result = container.add_chromosome("Nonexistent Genome", "1".to_string(), "ATCG".to_string());
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Genome 'Nonexistent Genome' not found");
    }

    #[test]
    fn test_add_duplicate_genome() {
        let mut container = GenomeContainer::new();
        let genome_metadata = GenomeDataLocal {
            description: Some("Test Genome".to_string()),
            ..Default::default()
        };
        container.add_genome("Test Genome".to_string(), genome_metadata.clone()).unwrap();
        let result = container.add_genome("Test Genome".to_string(), genome_metadata);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Genome 'Test Genome' already exists");
    }

    #[test]
    fn test_download_genome_files() -> Result<(), Box<dyn std::error::Error>> {
        let config_path = "assets/genomes.json";
        let mut container = GenomeContainer::from_json_file(config_path)?;

        container.download_genome_files("Saccharomyces cerevisiae S288c Ensembl 113")?;
        let metadata = container.get_genomes_metadata("Saccharomyces cerevisiae S288c Ensembl 113")?;
        assert!(metadata.downloaded);

        Ok(())
    }

    #[test]
    fn test_load_fasta_sequences() -> Result<(), Box<dyn std::error::Error>> {
        let config_path = "assets/genomes.json";
        let mut container = GenomeContainer::from_json_file(config_path)?;

        container.download_genome_files("Saccharomyces cerevisiae S288c Ensembl 113")?;
        container.load_fasta_sequences("Saccharomyces cerevisiae S288c Ensembl 113")?;

        let sequence = container.get_sequence("Saccharomyces cerevisiae S288c Ensembl 113", "I")?;
        assert!(sequence.len() > 10000);

        Ok(())
    }

    #[test]
    fn test_load_gff_annotations() -> Result<(), Box<dyn std::error::Error>> {
        let config_path = "assets/genomes.json";
        let mut container = GenomeContainer::from_json_file(config_path)?;

        container.download_genome_files("Saccharomyces cerevisiae S288c Ensembl 113")?;
        container.load_gff_annotations("Saccharomyces cerevisiae S288c Ensembl 113")?;

        Ok(())
    }

    #[test]
    fn test_gene_data() -> Result<(), Box<dyn std::error::Error>> {
        let config_path = "assets/genomes.json";
        let mut container = GenomeContainer::from_json_file(config_path)?;

        container.download_genome_files("Saccharomyces cerevisiae S288c Ensembl 113")?;
        container.load_gff_annotations("Saccharomyces cerevisiae S288c Ensembl 113")?;

        let gene = container.genes.get("YDR387C").ok_or("Gene 'YDR387C' not found")?;
        assert_eq!(gene.gene_name.as_deref(), Some("CIN10"));
        assert_eq!(gene.start, 1248154);
        assert_eq!(gene.end, 1249821);

        assert!(gene.exons.len() > 0);

        let exon_with_end_equal_to_gene_end = gene.exons.iter().any(|exon| exon.end == gene.end);
        assert!(exon_with_end_equal_to_gene_end);

        Ok(())
    }
}

