//! PrimerBank lookup records and the official web-search adapter.
//!
//! PrimerBank exposes a public HTML search form rather than a documented JSON
//! API. This module keeps that transport detail behind a typed report and also
//! accepts saved HTML, so tests and reproducible workflows do not require the
//! network. It deliberately retrieves individual search results and does not
//! mirror or bundle the PrimerBank database.

use scraper::{ElementRef, Html, Selector};
use serde::{Deserialize, Serialize};
use std::{fs, path::Path, time::Duration};

pub const PRIMERBANK_SEARCH_REPORT_SCHEMA: &str = "gentle.primerbank_search.v1";
pub const PRIMERBANK_CDNA_TEST_REPORT_SCHEMA: &str = "gentle.primerbank_cdna_test.v1";
pub const PRIMERBANK_HOME_URL: &str = "https://pga.mgh.harvard.edu/primerbank/";
pub const PRIMERBANK_USAGE_POLICY_URL: &str =
    "https://pga.mgh.harvard.edu/primerbank/citation.html";
pub const PRIMERBANK_SEARCH_URL: &str =
    "https://pga.mgh.harvard.edu/cgi-bin/primerbank/new_search2.cgi";
const PRIMERBANK_DETAIL_BASE_URL: &str = "https://pga.mgh.harvard.edu";
const PRIMERBANK_MAX_RESPONSE_BYTES: usize = 8 * 1024 * 1024;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PrimerBankQueryKind {
    GenbankAccession,
    NcbiProteinAccession,
    NcbiGeneId,
    PrimerbankId,
    #[default]
    NcbiGeneSymbol,
    Keyword,
}

impl PrimerBankQueryKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::GenbankAccession => "genbank_accession",
            Self::NcbiProteinAccession => "ncbi_protein_accession",
            Self::NcbiGeneId => "ncbi_gene_id",
            Self::PrimerbankId => "primerbank_id",
            Self::NcbiGeneSymbol => "ncbi_gene_symbol",
            Self::Keyword => "keyword",
        }
    }

    pub fn form_value(self) -> &'static str {
        match self {
            Self::GenbankAccession => "GenBank Accession",
            Self::NcbiProteinAccession => "NCBI Protein Accession",
            Self::NcbiGeneId => "NCBI Gene ID",
            Self::PrimerbankId => "PrimerBank ID",
            Self::NcbiGeneSymbol => "NCBI Gene Symbol",
            Self::Keyword => "Keyword",
        }
    }

    pub fn parse(raw: &str) -> Option<Self> {
        match normalize_token(raw).as_str() {
            "genbank" | "genbankaccession" => Some(Self::GenbankAccession),
            "protein" | "ncbiprotein" | "ncbiproteinaccession" => Some(Self::NcbiProteinAccession),
            "geneid" | "ncbigeneid" => Some(Self::NcbiGeneId),
            "primerbank" | "primerbankid" | "id" => Some(Self::PrimerbankId),
            "gene" | "symbol" | "genesymbol" | "ncbigenesymbol" => Some(Self::NcbiGeneSymbol),
            "keyword" | "description" => Some(Self::Keyword),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PrimerBankSpecies {
    All,
    #[default]
    Human,
    Mouse,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PrimerBankSpeciesMatchStatus {
    #[default]
    Unresolved,
    NotRequested,
    Matched,
    Mismatch,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct PrimerBankSpeciesCheck {
    pub requested_species: PrimerBankSpecies,
    pub status: PrimerBankSpeciesMatchStatus,
    pub matched_gene_count: usize,
    pub mismatched_gene_count: usize,
    pub unresolved_gene_count: usize,
    pub observed_species: Vec<String>,
}

impl PrimerBankSpecies {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::Human => "human",
            Self::Mouse => "mouse",
        }
    }

    pub fn form_value(self) -> &'static str {
        match self {
            Self::All => "All Species",
            Self::Human => "Human",
            Self::Mouse => "Mouse",
        }
    }

    pub fn parse(raw: &str) -> Option<Self> {
        match normalize_token(raw).as_str() {
            "all" | "allspecies" => Some(Self::All),
            "human" | "homosapiens" => Some(Self::Human),
            "mouse" | "musmusculus" => Some(Self::Mouse),
            _ => None,
        }
    }

    /// Compare this requested species with an organism label returned by a
    /// catalog or carried by an annotated target sequence.
    pub fn match_observed_label(self, observed: Option<&str>) -> PrimerBankSpeciesMatchStatus {
        species_match_status(self, observed)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct PrimerBankSearchRequest {
    pub query: String,
    pub query_kind: PrimerBankQueryKind,
    pub species: PrimerBankSpecies,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct PrimerBankPrimerRecord {
    pub role: String,
    pub sequence_5_to_3: String,
    pub length_nt: usize,
    pub tm_c: Option<f64>,
    pub location_raw: String,
    pub location_start_1based: Option<usize>,
    pub location_end_1based: Option<usize>,
    pub interval_start_1based: Option<usize>,
    pub interval_end_1based: Option<usize>,
    pub orientation: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct PrimerBankPrimerPairRecord {
    pub primerbank_id: String,
    pub amplicon_size_bp: Option<usize>,
    pub forward: PrimerBankPrimerRecord,
    pub reverse: PrimerBankPrimerRecord,
    pub detail_url: String,
    /// GENtle does not infer experimental validation from catalog presence.
    pub validation_status: String,
    pub coordinate_system: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct PrimerBankGeneRecord {
    pub ncbi_gene_id: Option<String>,
    pub genbank_accession: Option<String>,
    pub ncbi_protein_accession: Option<String>,
    pub species: Option<String>,
    pub species_match_status: PrimerBankSpeciesMatchStatus,
    pub coding_dna_length_bp: Option<usize>,
    pub gene_description: Option<String>,
    pub primer_pairs: Vec<PrimerBankPrimerPairRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct PrimerBankSearchReport {
    pub schema: String,
    pub query: PrimerBankSearchRequest,
    pub source_url: String,
    pub source_kind: String,
    pub usage_policy_url: String,
    pub gene_count: usize,
    pub primer_pair_count: usize,
    pub species_check: PrimerBankSpeciesCheck,
    pub genes: Vec<PrimerBankGeneRecord>,
    pub citations: Vec<String>,
    pub warnings: Vec<String>,
}

/// A catalog pair tested against GENtle's current transcript-derived cDNA model.
///
/// The nested cDNA result remains distinct from PrimerBank provenance: a
/// compatible amplicon does not establish whole-genome specificity or
/// experimental validation.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct PrimerBankCdnaTestReport {
    pub schema: String,
    pub primerbank_query: PrimerBankSearchRequest,
    pub primerbank_source_url: String,
    pub primerbank_usage_policy_url: String,
    pub primerbank_pair: PrimerBankPrimerPairRecord,
    pub expected_species: PrimerBankSpecies,
    pub primerbank_species: Option<String>,
    pub species_match_status: PrimerBankSpeciesMatchStatus,
    pub target_sequence_species: Option<String>,
    pub target_sequence_species_match_status: PrimerBankSpeciesMatchStatus,
    pub warnings: Vec<String>,
    pub interpretation: String,
    pub cdna_test: serde_json::Value,
}

impl PrimerBankSearchReport {
    pub fn primer_pairs(&self) -> impl Iterator<Item = &PrimerBankPrimerPairRecord> {
        self.genes.iter().flat_map(|gene| gene.primer_pairs.iter())
    }

    pub fn pair_by_id(&self, primerbank_id: &str) -> Option<&PrimerBankPrimerPairRecord> {
        self.primer_pairs().find(|pair| {
            pair.primerbank_id
                .eq_ignore_ascii_case(primerbank_id.trim())
        })
    }

    pub fn gene_and_pair_by_id(
        &self,
        primerbank_id: &str,
    ) -> Option<(&PrimerBankGeneRecord, &PrimerBankPrimerPairRecord)> {
        self.genes.iter().find_map(|gene| {
            gene.primer_pairs
                .iter()
                .find(|pair| {
                    pair.primerbank_id
                        .eq_ignore_ascii_case(primerbank_id.trim())
                })
                .map(|pair| (gene, pair))
        })
    }
}

#[derive(Debug)]
struct ParsedCell {
    text: String,
    hrefs: Vec<String>,
}

pub fn search_primerbank(
    request: &PrimerBankSearchRequest,
    source_html_path: Option<&Path>,
) -> Result<PrimerBankSearchReport, String> {
    validate_request(request)?;
    if let Some(path) = source_html_path {
        let html = fs::read_to_string(path).map_err(|error| {
            format!(
                "Could not read saved PrimerBank HTML '{}': {error}",
                path.display()
            )
        })?;
        return parse_primerbank_search_html(
            request,
            &format!("file://{}", path.display()),
            "saved_html",
            &html,
        );
    }

    let client = reqwest::blocking::Client::builder()
        .connect_timeout(Duration::from_secs(10))
        .timeout(Duration::from_secs(30))
        .user_agent(format!(
            "GENtle/{} PrimerBank lookup",
            env!("CARGO_PKG_VERSION")
        ))
        .build()
        .map_err(|error| format!("Could not build PrimerBank HTTP client: {error}"))?;
    let response = client
        .post(PRIMERBANK_SEARCH_URL)
        .form(&[
            ("selectBox", request.query_kind.form_value()),
            ("species", request.species.form_value()),
            ("searchBox", request.query.trim()),
            ("Submit", "Submit"),
        ])
        .send()
        .map_err(|error| format!("PrimerBank search request failed: {error}"))?;
    let status = response.status();
    if !status.is_success() {
        return Err(format!(
            "PrimerBank search returned HTTP {status} from '{PRIMERBANK_SEARCH_URL}'"
        ));
    }
    if response
        .content_length()
        .is_some_and(|length| length > PRIMERBANK_MAX_RESPONSE_BYTES as u64)
    {
        return Err(format!(
            "PrimerBank response exceeds the {} byte safety limit",
            PRIMERBANK_MAX_RESPONSE_BYTES
        ));
    }
    let bytes = response
        .bytes()
        .map_err(|error| format!("Could not read PrimerBank response: {error}"))?;
    if bytes.len() > PRIMERBANK_MAX_RESPONSE_BYTES {
        return Err(format!(
            "PrimerBank response exceeds the {} byte safety limit",
            PRIMERBANK_MAX_RESPONSE_BYTES
        ));
    }
    let html = String::from_utf8_lossy(&bytes);
    parse_primerbank_search_html(
        request,
        PRIMERBANK_SEARCH_URL,
        "primerbank_live_search",
        &html,
    )
}

pub fn parse_primerbank_search_html(
    request: &PrimerBankSearchRequest,
    source_url: &str,
    source_kind: &str,
    html: &str,
) -> Result<PrimerBankSearchReport, String> {
    validate_request(request)?;
    let document = Html::parse_document(html);
    let table_selector = selector("table")?;
    let row_selector = selector("tr")?;
    let cell_selector = selector("td, th")?;
    let link_selector = selector("a[href]")?;
    let document_text = normalized_text(document.root_element());
    if !document_text.to_ascii_lowercase().contains("primerbank") {
        return Err("Saved HTML does not look like a PrimerBank response".to_string());
    }

    let mut genes = Vec::<PrimerBankGeneRecord>::new();
    let mut current_gene_index: Option<usize> = None;
    for table in document.select(&table_selector) {
        // PrimerBank uses layout tables around the record tables. Only leaf
        // tables carry one semantic record; parsing an ancestor would count
        // the same gene twice and absorb all nested text into one row set.
        if table.select(&table_selector).next().is_some() {
            continue;
        }
        let rows = parse_table_rows(table, &row_selector, &cell_selector, &link_selector);
        let flattened = rows
            .iter()
            .flat_map(|row| row.iter().map(|cell| cell.text.as_str()))
            .collect::<Vec<_>>()
            .join(" ");
        if flattened.contains("Gene Description")
            && (flattened.contains("NCBI GeneID")
                || flattened.contains("GenBank Accession")
                || flattened.contains("Species"))
        {
            let gene = parse_gene_table(&rows);
            genes.push(gene);
            current_gene_index = Some(genes.len() - 1);
        } else if flattened.contains("PrimerBank ID")
            && flattened.contains("Forward Primer")
            && flattened.contains("Reverse Primer")
        {
            let pair = parse_pair_table(&rows)?;
            let index = if let Some(index) = current_gene_index {
                index
            } else {
                genes.push(PrimerBankGeneRecord::default());
                current_gene_index = Some(genes.len() - 1);
                genes.len() - 1
            };
            genes[index].primer_pairs.push(pair);
        }
    }
    genes.retain(|gene| {
        !gene.primer_pairs.is_empty()
            || gene.ncbi_gene_id.is_some()
            || gene.genbank_accession.is_some()
            || gene.gene_description.is_some()
    });
    for gene in &mut genes {
        gene.species_match_status = species_match_status(request.species, gene.species.as_deref());
    }
    let species_check = summarize_species_check(request.species, &genes);
    let primer_pair_count = genes.iter().map(|gene| gene.primer_pairs.len()).sum();
    let mut warnings = vec![
        "PrimerBank records are external catalog candidates; GENtle has not established compatibility with the current transcript annotation or specificity database."
            .to_string(),
        "Primer locations use PrimerBank's coding-sequence record (1-based, inclusive), not genomic or current-project coordinates."
            .to_string(),
        "Catalog presence is not reported as direct experimental validation; inspect the PrimerBank detail page and test the pair in the intended biological context."
            .to_string(),
        "Respect PrimerBank citation and usage policy; GENtle retrieves individual records and does not mirror the database."
            .to_string(),
    ];
    let lower_document = document_text.to_ascii_lowercase();
    if primer_pair_count == 0 {
        if lower_document.contains("no primer pair is found") {
            warnings.push("PrimerBank returned no matching primer pair.".to_string());
        } else {
            return Err(
                "PrimerBank response contained no recognizable primer-pair table; the upstream HTML layout may have changed"
                    .to_string(),
            );
        }
    }
    match species_check.status {
        PrimerBankSpeciesMatchStatus::Mismatch => warnings.push(format!(
            "Requested PrimerBank species '{}' does not match {} returned gene record(s); observed species: {}.",
            request.species.as_str(),
            species_check.mismatched_gene_count,
            species_check.observed_species.join(", ")
        )),
        PrimerBankSpeciesMatchStatus::Unresolved => warnings.push(format!(
            "Could not verify requested PrimerBank species '{}' for {} returned gene record(s).",
            request.species.as_str(),
            species_check.unresolved_gene_count
        )),
        PrimerBankSpeciesMatchStatus::NotRequested if primer_pair_count > 0 => warnings.push(
            "Species matching was not requested; use species=human or species=mouse before biological use."
                .to_string(),
        ),
        PrimerBankSpeciesMatchStatus::Matched | PrimerBankSpeciesMatchStatus::NotRequested => {}
    }
    Ok(PrimerBankSearchReport {
        schema: PRIMERBANK_SEARCH_REPORT_SCHEMA.to_string(),
        query: request.clone(),
        source_url: source_url.to_string(),
        source_kind: source_kind.to_string(),
        usage_policy_url: PRIMERBANK_USAGE_POLICY_URL.to_string(),
        gene_count: genes.len(),
        primer_pair_count,
        species_check,
        genes,
        citations: vec![
            "Spandidos A, Wang X, Wang H, Seed B. PrimerBank: a resource of human and mouse PCR primer pairs for gene expression detection and quantification. Nucleic Acids Res. 2010;38:D792-D799."
                .to_string(),
            "Spandidos A, Wang X, Wang H, Dragnev S, Thurber T, Seed B. A comprehensive collection of experimentally validated primers for Polymerase Chain Reaction quantitation of murine transcript abundance. BMC Genomics. 2008;9:633."
                .to_string(),
        ],
        warnings,
    })
}

fn validate_request(request: &PrimerBankSearchRequest) -> Result<(), String> {
    let query = request.query.trim();
    if query.is_empty() {
        return Err("PrimerBank query must not be empty".to_string());
    }
    if query.chars().count() > 500 {
        return Err("PrimerBank query must not exceed 500 characters".to_string());
    }
    if query.chars().any(char::is_control) {
        return Err("PrimerBank query must not contain control characters".to_string());
    }
    Ok(())
}

fn parse_gene_table(rows: &[Vec<ParsedCell>]) -> PrimerBankGeneRecord {
    PrimerBankGeneRecord {
        ncbi_gene_id: row_value(rows, "NCBI GeneID"),
        genbank_accession: row_value(rows, "GenBank Accession"),
        ncbi_protein_accession: row_value(rows, "NCBI Protein Accession"),
        species: row_value(rows, "Species"),
        species_match_status: PrimerBankSpeciesMatchStatus::Unresolved,
        coding_dna_length_bp: row_value(rows, "Coding DNA Length")
            .and_then(|value| value.parse::<usize>().ok()),
        gene_description: row_value(rows, "Gene Description"),
        primer_pairs: vec![],
    }
}

fn species_match_status(
    requested: PrimerBankSpecies,
    observed: Option<&str>,
) -> PrimerBankSpeciesMatchStatus {
    if requested == PrimerBankSpecies::All {
        return PrimerBankSpeciesMatchStatus::NotRequested;
    }
    let Some(observed) = observed.and_then(parse_species_label) else {
        return PrimerBankSpeciesMatchStatus::Unresolved;
    };
    if requested == observed {
        PrimerBankSpeciesMatchStatus::Matched
    } else {
        PrimerBankSpeciesMatchStatus::Mismatch
    }
}

fn parse_species_label(raw: &str) -> Option<PrimerBankSpecies> {
    match normalize_token(raw).as_str() {
        "human" | "homosapiens" => Some(PrimerBankSpecies::Human),
        "mouse" | "musmusculus" => Some(PrimerBankSpecies::Mouse),
        _ => None,
    }
}

fn summarize_species_check(
    requested: PrimerBankSpecies,
    genes: &[PrimerBankGeneRecord],
) -> PrimerBankSpeciesCheck {
    let mut observed_species = genes
        .iter()
        .filter_map(|gene| gene.species.clone())
        .collect::<Vec<_>>();
    observed_species.sort();
    observed_species.dedup();
    let matched_gene_count = genes
        .iter()
        .filter(|gene| gene.species_match_status == PrimerBankSpeciesMatchStatus::Matched)
        .count();
    let mismatched_gene_count = genes
        .iter()
        .filter(|gene| gene.species_match_status == PrimerBankSpeciesMatchStatus::Mismatch)
        .count();
    let unresolved_gene_count = genes
        .iter()
        .filter(|gene| gene.species_match_status == PrimerBankSpeciesMatchStatus::Unresolved)
        .count();
    let status = if requested == PrimerBankSpecies::All {
        PrimerBankSpeciesMatchStatus::NotRequested
    } else if mismatched_gene_count > 0 {
        PrimerBankSpeciesMatchStatus::Mismatch
    } else if genes.is_empty() || unresolved_gene_count > 0 {
        PrimerBankSpeciesMatchStatus::Unresolved
    } else {
        PrimerBankSpeciesMatchStatus::Matched
    };
    PrimerBankSpeciesCheck {
        requested_species: requested,
        status,
        matched_gene_count,
        mismatched_gene_count,
        unresolved_gene_count,
        observed_species,
    }
}

fn parse_pair_table(rows: &[Vec<ParsedCell>]) -> Result<PrimerBankPrimerPairRecord, String> {
    let primerbank_id = row_value(rows, "PrimerBank ID")
        .filter(|value| !value.is_empty())
        .ok_or_else(|| "PrimerBank primer-pair table is missing PrimerBank ID".to_string())?;
    let amplicon_size_bp =
        row_value(rows, "Amplicon Size").and_then(|value| value.parse::<usize>().ok());
    let forward = parse_primer_row(rows, "Forward Primer", "forward")?;
    let reverse = parse_primer_row(rows, "Reverse Primer", "reverse")?;
    let detail_href = rows
        .iter()
        .flat_map(|row| row.iter())
        .flat_map(|cell| cell.hrefs.iter())
        .find(|href| href.contains("displayDetail") && href.contains("primerID="))
        .cloned();
    let detail_url = detail_href
        .as_deref()
        .map(absolute_primerbank_url)
        .unwrap_or_else(|| {
            format!(
                "{PRIMERBANK_DETAIL_BASE_URL}/cgi-bin/primerbank/new_displayDetail2.cgi?primerID={primerbank_id}"
            )
        });
    Ok(PrimerBankPrimerPairRecord {
        primerbank_id,
        amplicon_size_bp,
        forward,
        reverse,
        detail_url,
        validation_status: "not_assessed_by_gentle".to_string(),
        coordinate_system: "primerbank_coding_sequence_1based_inclusive".to_string(),
    })
}

fn parse_primer_row(
    rows: &[Vec<ParsedCell>],
    label: &str,
    role: &str,
) -> Result<PrimerBankPrimerRecord, String> {
    let row = rows
        .iter()
        .find(|row| row.first().is_some_and(|cell| cell.text == label))
        .ok_or_else(|| format!("PrimerBank primer-pair table is missing {label}"))?;
    if row.len() < 5 {
        return Err(format!(
            "PrimerBank {label} row has {} cells; expected at least 5",
            row.len()
        ));
    }
    let sequence = normalize_dna(&row[1].text)?;
    let length_nt = row[2].text.parse::<usize>().unwrap_or(sequence.len());
    let tm_c = row[3].text.parse::<f64>().ok();
    let location_raw = row[4].text.clone();
    let (location_start_1based, location_end_1based) = parse_location(&location_raw);
    let (interval_start_1based, interval_end_1based) =
        match (location_start_1based, location_end_1based) {
            (Some(start), Some(end)) => (Some(start.min(end)), Some(start.max(end))),
            _ => (None, None),
        };
    Ok(PrimerBankPrimerRecord {
        role: role.to_string(),
        sequence_5_to_3: sequence,
        length_nt,
        tm_c,
        location_raw,
        location_start_1based,
        location_end_1based,
        interval_start_1based,
        interval_end_1based,
        orientation: if role == "reverse" {
            "reverse".to_string()
        } else {
            "forward".to_string()
        },
    })
}

fn parse_table_rows(
    table: ElementRef<'_>,
    row_selector: &Selector,
    cell_selector: &Selector,
    link_selector: &Selector,
) -> Vec<Vec<ParsedCell>> {
    table
        .select(row_selector)
        .map(|row| {
            row.select(cell_selector)
                .map(|cell| ParsedCell {
                    text: normalized_text(cell),
                    hrefs: cell
                        .select(link_selector)
                        .filter_map(|link| link.value().attr("href"))
                        .map(ToString::to_string)
                        .collect(),
                })
                .collect::<Vec<_>>()
        })
        .filter(|row| !row.is_empty())
        .collect()
}

fn row_value(rows: &[Vec<ParsedCell>], label: &str) -> Option<String> {
    rows.iter()
        .find(|row| row.first().is_some_and(|cell| cell.text == label))
        .and_then(|row| row.get(1))
        .map(|cell| cell.text.trim().to_string())
        .filter(|value| !value.is_empty())
}

fn normalized_text(element: ElementRef<'_>) -> String {
    element
        .text()
        .flat_map(str::split_whitespace)
        .collect::<Vec<_>>()
        .join(" ")
}

fn selector(raw: &str) -> Result<Selector, String> {
    Selector::parse(raw).map_err(|error| format!("Invalid internal HTML selector '{raw}': {error}"))
}

fn normalize_token(raw: &str) -> String {
    raw.chars()
        .filter(|ch| ch.is_ascii_alphanumeric())
        .flat_map(char::to_lowercase)
        .collect()
}

fn normalize_dna(raw: &str) -> Result<String, String> {
    let sequence = raw
        .chars()
        .filter(|ch| !ch.is_whitespace())
        .flat_map(char::to_uppercase)
        .collect::<String>();
    if sequence.is_empty()
        || !sequence
            .chars()
            .all(|ch| matches!(ch, 'A' | 'C' | 'G' | 'T'))
    {
        return Err(format!(
            "PrimerBank returned an invalid DNA primer sequence '{raw}'"
        ));
    }
    Ok(sequence)
}

fn parse_location(raw: &str) -> (Option<usize>, Option<usize>) {
    let mut parts = raw.splitn(2, '-');
    let start = parts.next().and_then(|value| value.trim().parse().ok());
    let end = parts.next().and_then(|value| value.trim().parse().ok());
    (start, end)
}

fn absolute_primerbank_url(href: &str) -> String {
    if href.starts_with("http://") || href.starts_with("https://") {
        href.to_string()
    } else if href.starts_with('/') {
        format!("{PRIMERBANK_DETAIL_BASE_URL}{href}")
    } else {
        format!("{PRIMERBANK_HOME_URL}{href}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_html() -> String {
        fs::read_to_string("test_files/fixtures/primerbank/synthetic_primerbank_search_result.html")
            .expect("synthetic PrimerBank fixture")
    }

    #[test]
    fn parses_synthetic_search_result_without_validation_claims() {
        let request = PrimerBankSearchRequest {
            query: "TOY1".to_string(),
            query_kind: PrimerBankQueryKind::NcbiGeneSymbol,
            species: PrimerBankSpecies::Human,
        };
        let report = parse_primerbank_search_html(
            &request,
            "synthetic://primerbank/search",
            "synthetic_fixture",
            &fixture_html(),
        )
        .expect("parse synthetic PrimerBank result");
        assert_eq!(report.schema, PRIMERBANK_SEARCH_REPORT_SCHEMA);
        assert_eq!(report.gene_count, 1);
        assert_eq!(report.primer_pair_count, 2);
        assert_eq!(
            report.species_check.status,
            PrimerBankSpeciesMatchStatus::Matched
        );
        assert_eq!(report.species_check.observed_species, vec!["Human"]);
        let gene = &report.genes[0];
        assert_eq!(gene.ncbi_gene_id.as_deref(), Some("12345"));
        assert_eq!(gene.genbank_accession.as_deref(), Some("NM_000001"));
        assert_eq!(
            gene.species_match_status,
            PrimerBankSpeciesMatchStatus::Matched
        );
        let pair = &gene.primer_pairs[0];
        assert_eq!(pair.primerbank_id, "100000001a1");
        assert_eq!(pair.amplicon_size_bp, Some(121));
        assert_eq!(pair.forward.sequence_5_to_3, "ACGTACGTACGTACGTACGTA");
        assert_eq!(pair.forward.location_start_1based, Some(101));
        assert_eq!(pair.forward.location_end_1based, Some(121));
        assert_eq!(pair.reverse.location_start_1based, Some(221));
        assert_eq!(pair.reverse.location_end_1based, Some(202));
        assert_eq!(pair.reverse.interval_start_1based, Some(202));
        assert_eq!(pair.reverse.interval_end_1based, Some(221));
        assert_eq!(pair.validation_status, "not_assessed_by_gentle");
        assert!(pair.detail_url.ends_with("primerID=100000001a1"));
    }

    #[test]
    fn reports_species_mismatch_without_hiding_catalog_rows() {
        let request = PrimerBankSearchRequest {
            query: "TOY1".to_string(),
            query_kind: PrimerBankQueryKind::NcbiGeneSymbol,
            species: PrimerBankSpecies::Mouse,
        };
        let report = parse_primerbank_search_html(
            &request,
            "synthetic://primerbank/search",
            "synthetic_fixture",
            &fixture_html(),
        )
        .expect("parse mismatching synthetic PrimerBank result");
        assert_eq!(report.primer_pair_count, 2);
        assert_eq!(
            report.species_check.status,
            PrimerBankSpeciesMatchStatus::Mismatch
        );
        assert_eq!(report.species_check.mismatched_gene_count, 1);
        assert_eq!(
            report.genes[0].species_match_status,
            PrimerBankSpeciesMatchStatus::Mismatch
        );
        assert!(
            report
                .warnings
                .iter()
                .any(|warning| warning.contains("does not match"))
        );
    }

    #[test]
    fn exact_id_shape_without_ncbi_gene_id_still_cross_checks_species() {
        let request = PrimerBankSearchRequest {
            query: "100000001a1".to_string(),
            query_kind: PrimerBankQueryKind::PrimerbankId,
            species: PrimerBankSpecies::Human,
        };
        let html = fixture_html().replace("NCBI GeneID", "Legacy Gene Identifier");
        let report = parse_primerbank_search_html(
            &request,
            "synthetic://primerbank/exact-id",
            "synthetic_fixture",
            &html,
        )
        .expect("parse exact-id response without an NCBI GeneID row");
        assert_eq!(report.genes[0].species.as_deref(), Some("Human"));
        assert_eq!(
            report.species_check.status,
            PrimerBankSpeciesMatchStatus::Matched
        );
    }

    #[test]
    fn parses_valid_empty_search_result() {
        let request = PrimerBankSearchRequest {
            query: "NO_MATCH".to_string(),
            ..PrimerBankSearchRequest::default()
        };
        let html = "<html><head><title>PrimerBank Search Result</title></head><body>No primer pair is found for NCBI Gene Symbol (Human): NO_MATCH</body></html>";
        let report = parse_primerbank_search_html(
            &request,
            "synthetic://primerbank/empty",
            "synthetic_fixture",
            html,
        )
        .expect("valid empty report");
        assert_eq!(report.primer_pair_count, 0);
        assert!(
            report
                .warnings
                .iter()
                .any(|warning| warning.contains("no matching"))
        );
    }

    #[test]
    fn query_kind_and_species_accept_cli_aliases() {
        assert_eq!(
            PrimerBankQueryKind::parse("gene-symbol"),
            Some(PrimerBankQueryKind::NcbiGeneSymbol)
        );
        assert_eq!(
            PrimerBankQueryKind::parse("PrimerBank-ID"),
            Some(PrimerBankQueryKind::PrimerbankId)
        );
        assert_eq!(
            PrimerBankSpecies::parse("Homo sapiens"),
            Some(PrimerBankSpecies::Human)
        );
    }
}
