//! Optional Ensembl Regulation source, snapshot, index, and overlap helpers.
//!
//! Provider API rows are normalized into a deterministic local TSV before a
//! sparse byte-offset index is built. Neither file is opened during ordinary
//! startup; callers opt into resource installation and querying explicitly.

use crate::digest_utils::sha256_file_hex;
use gentle_protocol::{
    ENSEMBL_REGULATION_INTERVAL_INDEX_SCHEMA, ENSEMBL_REGULATION_SOURCE_CATALOG_SCHEMA,
    EnsemblRegulationBinCheckpoint, EnsemblRegulationChromosomeIndex, EnsemblRegulationInterval,
    EnsemblRegulationIntervalIndex, EnsemblRegulationSourceCatalog,
    EnsemblRegulationSourceDescriptor,
};
use serde::Deserialize;
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
};

pub const ENSEMBL_REGULATION_INTERVALS_SCHEMA: &str = "gentle.ensembl_regulation_intervals.v1";
pub const DEFAULT_ENSEMBL_REGULATION_BIN_SIZE_BP: u64 = 65_536;
pub const ENSEMBL_REGULATION_2026_08_GRCH38_SOURCE_ID: &str = "ensembl_regulation_2026_08_grch38";
pub const ENSEMBL_REGULATION_2026_08_GRCM39_SOURCE_ID: &str = "ensembl_regulation_2026_08_grcm39";
pub const ENSEMBL_REGULATION_DATA_ACCESS_URL: &str = "https://regulation.ensembl.org/data_access";
pub const ENSEMBL_REGULATION_PRIMARY_ANALYSIS_URL: &str =
    "https://regulation.ensembl.org/epigenomes";

fn source(
    source_id: &str,
    species: &str,
    taxon_id: u32,
    assembly_name: &str,
    assembly_aliases: &[&str],
    assembly_accession: &str,
    gene_annotation_release: &str,
    activity_url: &str,
    browser_species_slug: &str,
) -> EnsemblRegulationSourceDescriptor {
    let annotation_release = "2026-08";
    let api_version = "v0.15";
    EnsemblRegulationSourceDescriptor {
        source_id: source_id.to_string(),
        provider: "Ensembl Regulation".to_string(),
        annotation_release: annotation_release.to_string(),
        annotation_api_version: api_version.to_string(),
        pipeline_version: "2.1".to_string(),
        gene_annotation_source: "GENCODE basic".to_string(),
        gene_annotation_release: gene_annotation_release.to_string(),
        species_scientific_name: species.to_string(),
        taxon_id,
        assembly_name: assembly_name.to_string(),
        assembly_aliases: assembly_aliases
            .iter()
            .map(|value| (*value).to_string())
            .collect(),
        assembly_accession: assembly_accession.to_string(),
        feature_types: vec![
            "promoter".to_string(),
            "enhancer".to_string(),
            "open_chromatin_region".to_string(),
            "ctcf".to_string(),
            "emar".to_string(),
        ],
        annotation_api_url: format!(
            "https://regulation.ensembl.org/api/annotation/{api_version}/regulatory-features/ensembl/{annotation_release}/assembly/{assembly_name}"
        ),
        regulatory_activity_url: Some(activity_url.to_string()),
        data_access_url: ENSEMBL_REGULATION_DATA_ACCESS_URL.to_string(),
        primary_analysis_url: ENSEMBL_REGULATION_PRIMARY_ANALYSIS_URL.to_string(),
        browser_species_slug: browser_species_slug.to_string(),
        feature_page_url_template: format!(
            "https://regulation.ensembl.org/regulatory_features/{browser_species_slug}/{{FEATURE_ID}}"
        ),
        coordinate_system:
            "Ensembl API 1-based inclusive; GENtle normalized TSV 0-based half-open".to_string(),
        scope_note: "Release-bound Ensembl regulatory annotation. Associated genes are provider annotations, not proof that a feature regulates a named gene."
            .to_string(),
        activity_note: "The separate regulatory-activity matrix is not imported by the annotation snapshot route and must not be inferred from feature overlap."
            .to_string(),
        primary_signal_note: "Epigenome-specific peak and signal files are separate optional evidence. Quantitative signal must retain assay, biosample, replicate, units, and source provenance."
            .to_string(),
    }
}

pub fn source_catalog() -> EnsemblRegulationSourceCatalog {
    EnsemblRegulationSourceCatalog {
        schema: ENSEMBL_REGULATION_SOURCE_CATALOG_SCHEMA.to_string(),
        sources: vec![
            source(
                ENSEMBL_REGULATION_2026_08_GRCH38_SOURCE_ID,
                "Homo sapiens",
                9606,
                "GRCh38",
                &["GRCh38", "hg38"],
                "GCA_000001405.29",
                "50",
                "https://regulation.ensembl.org/api/annotation/v0.15/files/download/2026-08/homo_sapiens/GRCh38/Homo_sapiens.GRCh38.regulatory_activity.tsv.gz",
                "homo_sapiens",
            ),
            source(
                ENSEMBL_REGULATION_2026_08_GRCM39_SOURCE_ID,
                "Mus musculus",
                10090,
                "GRCm39",
                &["GRCm39", "mm39"],
                "GCA_000001635.9",
                "M39",
                "https://regulation.ensembl.org/api/annotation/v0.15/files/download/2026-08/mus_musculus/GRCm39/Mus_musculus.GRCm39.regulatory_activity.tsv.gz",
                "mus_musculus",
            ),
        ],
        notes: vec![
            "No Ensembl Regulation resource is fetched, opened, or indexed until an explicit resource command is run."
                .to_string(),
            "Ensembl annotation, regulatory activity, and primary peak/signal files are separate evidence layers."
                .to_string(),
            "Ensembl mouse GRCm39/mm39 is not assembly-compatible with SCREEN Registry V4 mouse GRCm38/mm10."
                .to_string(),
            "The source catalog is pinned; GENtle does not silently replace it with the provider's latest release."
                .to_string(),
        ],
    }
}

pub fn source_descriptor(source_id: &str) -> Result<EnsemblRegulationSourceDescriptor, String> {
    let source_id = source_id.trim();
    source_catalog()
        .sources
        .into_iter()
        .find(|row| row.source_id == source_id)
        .ok_or_else(|| {
            format!(
                "Unknown Ensembl Regulation source '{source_id}'; available source ids: {}",
                source_catalog()
                    .sources
                    .iter()
                    .map(|row| row.source_id.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        })
}

/// Builds one canonical provider feature URL from the pinned source contract.
pub fn feature_page_url(
    source: &EnsemblRegulationSourceDescriptor,
    feature_id: &str,
) -> Result<String, String> {
    let catalog_source = source_descriptor(&source.source_id)?;
    if source.browser_species_slug != catalog_source.browser_species_slug
        || source.feature_page_url_template != catalog_source.feature_page_url_template
    {
        return Err(format!(
            "Ensembl Regulation source '{}' does not match the pinned feature-browser contract",
            source.source_id
        ));
    }
    let feature_id = feature_id.trim();
    let safe_feature_id = feature_id.starts_with("ENSR")
        && feature_id.len() <= 128
        && feature_id
            .chars()
            .all(|ch| ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.'));
    if !safe_feature_id {
        return Err(format!(
            "Ensembl Regulation feature id '{feature_id}' is not safe for a provider URL"
        ));
    }
    if catalog_source
        .feature_page_url_template
        .matches("{FEATURE_ID}")
        .count()
        != 1
    {
        return Err(format!(
            "Ensembl Regulation source '{}' has an invalid feature URL template",
            source.source_id
        ));
    }
    let url = catalog_source
        .feature_page_url_template
        .replace("{FEATURE_ID}", feature_id);
    if !url.starts_with("https://regulation.ensembl.org/") {
        return Err(format!(
            "Ensembl Regulation source '{}' produced an unsafe feature URL",
            source.source_id
        ));
    }
    Ok(url)
}

fn path_token(raw: &str) -> String {
    raw.chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() {
                ch.to_ascii_lowercase()
            } else {
                '.'
            }
        })
        .collect::<String>()
        .split('.')
        .filter(|part| !part.is_empty())
        .collect::<Vec<_>>()
        .join(".")
}

pub fn default_intervals_path(source: &EnsemblRegulationSourceDescriptor) -> String {
    format!(
        "data/resources/ensembl-regulation.{}.{}.intervals.tsv",
        path_token(&source.annotation_release),
        path_token(&source.assembly_name)
    )
}

pub fn default_index_path(source: &EnsemblRegulationSourceDescriptor) -> String {
    format!("{}.interval-index.json", default_intervals_path(source))
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub(crate) struct ApiResponse {
    pub release: String,
    pub regions: Vec<ApiRegionRows>,
    pub limit: usize,
    pub offset: Option<usize>,
    pub total_count: usize,
}

impl Default for ApiResponse {
    fn default() -> Self {
        Self {
            release: String::new(),
            regions: vec![],
            limit: 0,
            offset: None,
            total_count: 0,
        }
    }
}

#[derive(Debug, Clone, Deserialize, Default)]
#[serde(default)]
pub(crate) struct ApiRegionRows {
    pub region: ApiRegion,
    pub features: Vec<ApiFeature>,
}

#[derive(Debug, Clone, Deserialize, Default)]
#[serde(default)]
pub(crate) struct ApiRegion {
    pub name: String,
    pub length: u64,
    pub coordinate_system: String,
}

#[derive(Debug, Clone, Deserialize, Default)]
#[serde(default)]
pub(crate) struct ApiFeature {
    pub id: String,
    pub feature_type: String,
    pub start: u64,
    pub end: u64,
    pub strand: String,
    pub extended_start: Option<u64>,
    pub extended_end: Option<u64>,
    pub associated_genes: Vec<String>,
    pub associated_gene_names: Vec<String>,
}

fn normalized_strings(values: &[String]) -> Vec<String> {
    values
        .iter()
        .map(|value| value.trim())
        .filter(|value| !value.is_empty())
        .map(str::to_string)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn api_coordinate_start(value: Option<u64>, label: &str, id: &str) -> Result<Option<u64>, String> {
    value
        .map(|value| {
            value
                .checked_sub(1)
                .ok_or_else(|| format!("Ensembl Regulation feature '{id}' has invalid {label}=0"))
        })
        .transpose()
}

pub(crate) fn intervals_from_api_response(
    response: &ApiResponse,
    source: &EnsemblRegulationSourceDescriptor,
) -> Result<Vec<EnsemblRegulationInterval>, String> {
    if response.release != source.annotation_release {
        return Err(format!(
            "Ensembl Regulation API returned release '{}' for pinned source '{}' release '{}'; refusing a stale or redirected snapshot",
            response.release, source.source_id, source.annotation_release
        ));
    }
    let accepted_types = source
        .feature_types
        .iter()
        .map(|value| value.to_ascii_lowercase())
        .collect::<BTreeSet<_>>();
    let mut rows = Vec::new();
    for region_rows in &response.regions {
        let chromosome = region_rows.region.name.trim();
        if chromosome.is_empty() || region_rows.region.coordinate_system != "chromosome" {
            return Err(format!(
                "Ensembl Regulation API returned unsupported region '{}'/coordinate system '{}'",
                chromosome, region_rows.region.coordinate_system
            ));
        }
        for feature in &region_rows.features {
            if feature.id.trim().is_empty()
                || feature.start == 0
                || feature.end < feature.start
                || !accepted_types.contains(&feature.feature_type.to_ascii_lowercase())
            {
                return Err(format!(
                    "Ensembl Regulation API returned invalid feature '{}' ({}, {}..{})",
                    feature.id, feature.feature_type, feature.start, feature.end
                ));
            }
            if region_rows.region.length > 0 && feature.end > region_rows.region.length {
                return Err(format!(
                    "Ensembl Regulation feature '{}' ends beyond chromosome '{}' length {}",
                    feature.id, chromosome, region_rows.region.length
                ));
            }
            let strand = match feature.strand.as_str() {
                "forward" | "+" => Some('+'),
                "reverse" | "-" => Some('-'),
                "independent" | "." | "" => None,
                other => {
                    return Err(format!(
                        "Ensembl Regulation feature '{}' has unsupported strand '{}'",
                        feature.id, other
                    ));
                }
            };
            rows.push(EnsemblRegulationInterval {
                chromosome: chromosome.to_string(),
                start_0based: feature.start - 1,
                end_0based_exclusive: feature.end,
                feature_id: feature.id.trim().to_string(),
                feature_type: feature.feature_type.trim().to_string(),
                strand,
                extended_start_0based: api_coordinate_start(
                    feature.extended_start,
                    "extended_start",
                    &feature.id,
                )?,
                extended_end_0based_exclusive: feature.extended_end,
                associated_gene_ids: normalized_strings(&feature.associated_genes),
                associated_gene_names: normalized_strings(&feature.associated_gene_names),
                source_line_number: 0,
            });
        }
    }
    Ok(rows)
}

pub(crate) fn write_intervals_header<W: Write>(
    writer: &mut W,
    source: &EnsemblRegulationSourceDescriptor,
) -> Result<(), String> {
    writeln!(writer, "# {ENSEMBL_REGULATION_INTERVALS_SCHEMA}")
        .and_then(|_| writeln!(writer, "# source_id={}", source.source_id))
        .and_then(|_| {
            writeln!(
                writer,
                "# annotation_release={}",
                source.annotation_release
            )
        })
        .and_then(|_| writeln!(writer, "# assembly_name={}", source.assembly_name))
        .and_then(|_| {
            writeln!(
                writer,
                "# assembly_accession={}",
                source.assembly_accession
            )
        })
        .and_then(|_| {
            writeln!(
                writer,
                "# chromosome\tstart_0based\tend_0based_exclusive\tfeature_id\tfeature_type\tstrand\textended_start_0based\textended_end_0based_exclusive\tassociated_gene_ids_json\tassociated_gene_names_json"
            )
        })
        .map_err(|error| format!("Could not write Ensembl Regulation interval header: {error}"))
}

pub(crate) fn write_interval<W: Write>(
    writer: &mut W,
    row: &EnsemblRegulationInterval,
) -> Result<(), String> {
    let gene_ids = serde_json::to_string(&row.associated_gene_ids)
        .map_err(|error| format!("Could not serialize associated gene ids: {error}"))?;
    let gene_names = serde_json::to_string(&row.associated_gene_names)
        .map_err(|error| format!("Could not serialize associated gene names: {error}"))?;
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        row.chromosome,
        row.start_0based,
        row.end_0based_exclusive,
        row.feature_id,
        row.feature_type,
        row.strand
            .map(|value| value.to_string())
            .unwrap_or_else(|| ".".to_string()),
        row.extended_start_0based
            .map(|value| value.to_string())
            .unwrap_or_else(|| ".".to_string()),
        row.extended_end_0based_exclusive
            .map(|value| value.to_string())
            .unwrap_or_else(|| ".".to_string()),
        gene_ids,
        gene_names,
    )
    .map_err(|error| format!("Could not write Ensembl Regulation interval: {error}"))
}

fn parse_optional_u64(raw: &str, label: &str, line_number: u64) -> Result<Option<u64>, String> {
    if raw == "." || raw.is_empty() {
        return Ok(None);
    }
    raw.parse::<u64>().map(Some).map_err(|error| {
        format!("Invalid Ensembl Regulation {label} on line {line_number}: {error}")
    })
}

pub fn parse_interval_line(
    line: &str,
    source_line_number: u64,
) -> Result<Option<EnsemblRegulationInterval>, String> {
    let trimmed = line.trim_end_matches(['\r', '\n']);
    if trimmed.trim().is_empty() || trimmed.starts_with('#') {
        return Ok(None);
    }
    let fields = trimmed.split('\t').collect::<Vec<_>>();
    if fields.len() != 10 {
        return Err(format!(
            "Ensembl Regulation interval line {source_line_number} has {} fields; expected exactly 10",
            fields.len()
        ));
    }
    let start_0based = fields[1].parse::<u64>().map_err(|error| {
        format!("Invalid Ensembl Regulation start on line {source_line_number}: {error}")
    })?;
    let end_0based_exclusive = fields[2].parse::<u64>().map_err(|error| {
        format!("Invalid Ensembl Regulation end on line {source_line_number}: {error}")
    })?;
    if fields[0].trim().is_empty()
        || fields[3].trim().is_empty()
        || fields[4].trim().is_empty()
        || end_0based_exclusive <= start_0based
    {
        return Err(format!(
            "Invalid Ensembl Regulation interval on line {source_line_number}"
        ));
    }
    let strand = match fields[5] {
        "+" => Some('+'),
        "-" => Some('-'),
        "." | "" => None,
        other => {
            return Err(format!(
                "Invalid Ensembl Regulation strand '{other}' on line {source_line_number}"
            ));
        }
    };
    let associated_gene_ids = serde_json::from_str::<Vec<String>>(fields[8]).map_err(|error| {
        format!("Invalid associated gene-id JSON on line {source_line_number}: {error}")
    })?;
    let associated_gene_names =
        serde_json::from_str::<Vec<String>>(fields[9]).map_err(|error| {
            format!("Invalid associated gene-name JSON on line {source_line_number}: {error}")
        })?;
    Ok(Some(EnsemblRegulationInterval {
        chromosome: fields[0].trim().to_string(),
        start_0based,
        end_0based_exclusive,
        feature_id: fields[3].trim().to_string(),
        feature_type: fields[4].trim().to_string(),
        strand,
        extended_start_0based: parse_optional_u64(fields[6], "extended start", source_line_number)?,
        extended_end_0based_exclusive: parse_optional_u64(
            fields[7],
            "extended end",
            source_line_number,
        )?,
        associated_gene_ids: normalized_strings(&associated_gene_ids),
        associated_gene_names: normalized_strings(&associated_gene_names),
        source_line_number,
    }))
}

fn chromosome_aliases(chromosome: &str) -> Vec<String> {
    let mut aliases = BTreeSet::new();
    aliases.insert(chromosome.to_string());
    let without_chr = chromosome
        .strip_prefix("chr")
        .or_else(|| chromosome.strip_prefix("CHR"));
    if let Some(value) = without_chr {
        aliases.insert(value.to_string());
        if value.eq_ignore_ascii_case("M") {
            aliases.insert("MT".to_string());
            aliases.insert("chrMT".to_string());
        }
    } else {
        aliases.insert(format!("chr{chromosome}"));
    }
    aliases.into_iter().collect()
}

fn write_index(path: &str, index: &EnsemblRegulationIntervalIndex) -> Result<(), String> {
    let output = Path::new(path);
    if let Some(parent) = output
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent).map_err(|error| {
            format!(
                "Could not create Ensembl Regulation index directory '{}': {error}",
                parent.display()
            )
        })?;
    }
    let file = File::create(output)
        .map_err(|error| format!("Could not create Ensembl Regulation index '{path}': {error}"))?;
    let mut writer = BufWriter::new(file);
    serde_json::to_writer_pretty(&mut writer, index).map_err(|error| {
        format!("Could not serialize Ensembl Regulation index '{path}': {error}")
    })?;
    writer
        .write_all(b"\n")
        .and_then(|_| writer.flush())
        .map_err(|error| format!("Could not finish Ensembl Regulation index '{path}': {error}"))
}

pub fn prepare_interval_index(
    intervals_path: &str,
    output_path: &str,
    source: EnsemblRegulationSourceDescriptor,
) -> Result<EnsemblRegulationIntervalIndex, String> {
    let file = File::open(intervals_path).map_err(|error| {
        format!("Could not open Ensembl Regulation intervals '{intervals_path}': {error}")
    })?;
    let intervals_size_bytes = file
        .metadata()
        .map_err(|error| {
            format!("Could not inspect Ensembl Regulation intervals '{intervals_path}': {error}")
        })?
        .len();
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let mut byte_offset = 0u64;
    let mut source_line_number = 0u64;
    let mut row_count = 0u64;
    let mut feature_type_counts = BTreeMap::<String, u64>::new();
    let mut chromosome_bins = BTreeMap::<String, BTreeMap<u64, (u64, u64)>>::new();
    let mut chromosome_meta = BTreeMap::<String, (u64, u64, u64)>::new();
    let mut seen_chromosomes = BTreeSet::<String>::new();
    let mut current_chromosome: Option<String> = None;
    let mut previous_start = 0u64;
    let mut header_identity = BTreeMap::<String, String>::new();
    let accepted_types = source
        .feature_types
        .iter()
        .map(|value| value.to_ascii_lowercase())
        .collect::<BTreeSet<_>>();

    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line).map_err(|error| {
            format!("Could not read Ensembl Regulation intervals '{intervals_path}': {error}")
        })?;
        if bytes_read == 0 {
            break;
        }
        source_line_number = source_line_number.saturating_add(1);
        let next_byte_offset = byte_offset.saturating_add(bytes_read as u64);
        if let Some(raw) = line.trim().strip_prefix('#')
            && let Some((key, value)) = raw.trim().split_once('=')
        {
            header_identity.insert(key.trim().to_string(), value.trim().to_string());
        }
        let Some(interval) = parse_interval_line(&line, source_line_number)? else {
            byte_offset = next_byte_offset;
            continue;
        };
        if !accepted_types.contains(&interval.feature_type.to_ascii_lowercase()) {
            return Err(format!(
                "Ensembl Regulation interval line {} has type '{}' outside source '{}' declared types {}",
                source_line_number,
                interval.feature_type,
                source.source_id,
                source.feature_types.join(", ")
            ));
        }
        match current_chromosome.as_deref() {
            Some(current) if current == interval.chromosome => {
                if interval.start_0based < previous_start {
                    return Err(format!(
                        "Ensembl Regulation intervals are not start-sorted on chromosome '{}' at line {}",
                        interval.chromosome, source_line_number
                    ));
                }
            }
            Some(_) => {
                if seen_chromosomes.contains(&interval.chromosome) {
                    return Err(format!(
                        "Ensembl Regulation chromosome '{}' is split across non-contiguous blocks (line {})",
                        interval.chromosome, source_line_number
                    ));
                }
                seen_chromosomes.insert(interval.chromosome.clone());
                current_chromosome = Some(interval.chromosome.clone());
            }
            None => {
                seen_chromosomes.insert(interval.chromosome.clone());
                current_chromosome = Some(interval.chromosome.clone());
            }
        }
        previous_start = interval.start_0based;
        let first_bin = interval.start_0based / DEFAULT_ENSEMBL_REGULATION_BIN_SIZE_BP;
        let last_bin = interval.end_0based_exclusive.saturating_sub(1)
            / DEFAULT_ENSEMBL_REGULATION_BIN_SIZE_BP;
        let bins = chromosome_bins
            .entry(interval.chromosome.clone())
            .or_default();
        for bin_index in first_bin..=last_bin {
            bins.entry(bin_index)
                .and_modify(|entry| {
                    if byte_offset < entry.0 {
                        *entry = (byte_offset, source_line_number);
                    }
                })
                .or_insert((byte_offset, source_line_number));
        }
        chromosome_meta
            .entry(interval.chromosome.clone())
            .and_modify(|entry| {
                entry.1 = next_byte_offset;
                entry.2 = entry.2.saturating_add(1);
            })
            .or_insert((byte_offset, next_byte_offset, 1));
        *feature_type_counts
            .entry(interval.feature_type)
            .or_insert(0) += 1;
        row_count = row_count.saturating_add(1);
        byte_offset = next_byte_offset;
    }
    if row_count == 0 {
        return Err(format!(
            "Ensembl Regulation intervals '{intervals_path}' contain no data rows"
        ));
    }
    for (key, expected) in [
        ("source_id", source.source_id.as_str()),
        ("annotation_release", source.annotation_release.as_str()),
        ("assembly_name", source.assembly_name.as_str()),
        ("assembly_accession", source.assembly_accession.as_str()),
    ] {
        let observed = header_identity.get(key).map(String::as_str).unwrap_or("");
        if observed != expected {
            return Err(format!(
                "Ensembl Regulation intervals '{}' declare {key}='{}', expected '{}'; verify release/species/assembly before indexing",
                intervals_path, observed, expected
            ));
        }
    }
    let mut chromosome_alias_map = BTreeMap::new();
    let mut chromosomes = BTreeMap::new();
    for (chromosome, bins) in chromosome_bins {
        for alias in chromosome_aliases(&chromosome) {
            chromosome_alias_map.insert(alias, chromosome.clone());
        }
        let (first_row_byte, end_byte_exclusive, chromosome_row_count) = chromosome_meta
            .get(&chromosome)
            .copied()
            .unwrap_or_default();
        chromosomes.insert(
            chromosome.clone(),
            EnsemblRegulationChromosomeIndex {
                chromosome,
                row_count: chromosome_row_count,
                first_row_byte,
                end_byte_exclusive,
                bins: bins
                    .into_iter()
                    .map(|(bin_index, (scan_start_byte, scan_start_line_number))| {
                        EnsemblRegulationBinCheckpoint {
                            bin_index,
                            scan_start_byte,
                            scan_start_line_number,
                        }
                    })
                    .collect(),
            },
        );
    }
    let intervals_sha256 = format!(
        "sha256:{}",
        sha256_file_hex(Path::new(intervals_path)).map_err(|error| {
            format!("Could not hash Ensembl Regulation intervals '{intervals_path}': {error}")
        })?
    );
    let index = EnsemblRegulationIntervalIndex {
        schema: ENSEMBL_REGULATION_INTERVAL_INDEX_SCHEMA.to_string(),
        source,
        intervals_path_hint: intervals_path.to_string(),
        intervals_file_name: Path::new(intervals_path)
            .file_name()
            .and_then(|value| value.to_str())
            .unwrap_or_default()
            .to_string(),
        intervals_sha256,
        intervals_size_bytes,
        row_count,
        chromosome_count: chromosomes.len(),
        feature_type_counts,
        bin_size_bp: DEFAULT_ENSEMBL_REGULATION_BIN_SIZE_BP,
        chromosomes,
        chromosome_aliases: chromosome_alias_map,
        warnings: vec![],
    };
    write_index(output_path, &index)?;
    Ok(index)
}

pub fn read_interval_index(path: &str) -> Result<EnsemblRegulationIntervalIndex, String> {
    let text = fs::read_to_string(path)
        .map_err(|error| format!("Could not read Ensembl Regulation index '{path}': {error}"))?;
    let mut index = serde_json::from_str::<EnsemblRegulationIntervalIndex>(&text)
        .map_err(|error| format!("Invalid Ensembl Regulation index '{path}': {error}"))?;
    if index.schema != ENSEMBL_REGULATION_INTERVAL_INDEX_SCHEMA {
        return Err(format!(
            "Unsupported Ensembl Regulation index schema '{}' in '{}'; expected '{}'",
            index.schema, path, ENSEMBL_REGULATION_INTERVAL_INDEX_SCHEMA
        ));
    }
    if index.bin_size_bp == 0 || index.intervals_sha256.trim().is_empty() {
        return Err(format!(
            "Ensembl Regulation index '{path}' lacks a valid bin size or content digest"
        ));
    }
    let catalog_source = source_descriptor(&index.source.source_id).map_err(|error| {
        format!("Ensembl Regulation index '{path}' does not identify a supported source: {error}")
    })?;
    // Indexes created before the browser-link contract remain usable. Empty
    // serde-defaulted fields are hydrated only from the same pinned source;
    // non-empty mismatches still fail closed below.
    if index.source.browser_species_slug.trim().is_empty()
        && index.source.feature_page_url_template.trim().is_empty()
    {
        index.source.browser_species_slug = catalog_source.browser_species_slug.clone();
        index.source.feature_page_url_template = catalog_source.feature_page_url_template.clone();
    }
    if index.source != catalog_source {
        return Err(format!(
            "Ensembl Regulation index '{path}' source descriptor does not match GENtle's '{}' catalog entry; rebuild the index before use",
            index.source.source_id
        ));
    }
    Ok(index)
}

pub fn resolve_intervals_path(
    index_path: &str,
    index: &EnsemblRegulationIntervalIndex,
    override_path: Option<&str>,
) -> Result<PathBuf, String> {
    if let Some(path) = override_path
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        let path = PathBuf::from(path);
        if path.is_file() {
            return Ok(path);
        }
        return Err(format!(
            "Explicit Ensembl Regulation interval override '{}' does not exist",
            path.display()
        ));
    }
    if !index.intervals_file_name.is_empty() {
        let sibling = Path::new(index_path)
            .parent()
            .unwrap_or_else(|| Path::new(""))
            .join(&index.intervals_file_name);
        if sibling.is_file() {
            return Ok(sibling);
        }
    }
    let hint = PathBuf::from(&index.intervals_path_hint);
    if hint.is_file() {
        return Ok(hint);
    }
    Err(format!(
        "Could not locate Ensembl Regulation intervals for index '{}'; expected sibling '{}' or recorded path '{}'. Pass --intervals PATH after relocation.",
        index_path, index.intervals_file_name, index.intervals_path_hint
    ))
}

fn resolve_chromosome<'a>(
    index: &'a EnsemblRegulationIntervalIndex,
    requested: &str,
) -> Option<&'a str> {
    if let Some(canonical) = index.chromosome_aliases.get(requested) {
        return Some(canonical.as_str());
    }
    index
        .chromosome_aliases
        .iter()
        .find(|(alias, _)| alias.eq_ignore_ascii_case(requested))
        .map(|(_, canonical)| canonical.as_str())
}

pub fn validate_intervals_identity(
    path: &Path,
    index: &EnsemblRegulationIntervalIndex,
) -> Result<(), String> {
    let metadata = fs::metadata(path).map_err(|error| {
        format!(
            "Could not inspect Ensembl Regulation intervals '{}': {error}",
            path.display()
        )
    })?;
    if metadata.len() != index.intervals_size_bytes {
        return Err(format!(
            "Ensembl Regulation intervals '{}' have {} bytes, but the index binds {} bytes; rebuild the index",
            path.display(),
            metadata.len(),
            index.intervals_size_bytes
        ));
    }
    let observed = format!(
        "sha256:{}",
        sha256_file_hex(path).map_err(|error| {
            format!(
                "Could not hash Ensembl Regulation intervals '{}': {error}",
                path.display()
            )
        })?
    );
    if observed != index.intervals_sha256 {
        return Err(format!(
            "Ensembl Regulation intervals '{}' content digest '{}' does not match index digest '{}'; rebuild the index",
            path.display(),
            observed,
            index.intervals_sha256
        ));
    }
    Ok(())
}

pub fn query_overlaps(
    index: &EnsemblRegulationIntervalIndex,
    intervals_path: &Path,
    chromosome: &str,
    start_0based: u64,
    end_0based_exclusive: u64,
    requested_types: &[String],
) -> Result<Vec<EnsemblRegulationInterval>, String> {
    if end_0based_exclusive <= start_0based {
        return Ok(vec![]);
    }
    let Some(canonical) = resolve_chromosome(index, chromosome) else {
        return Ok(vec![]);
    };
    let Some(chromosome_index) = index.chromosomes.get(canonical) else {
        return Ok(vec![]);
    };
    let first_bin = start_0based / index.bin_size_bp;
    let last_bin = end_0based_exclusive.saturating_sub(1) / index.bin_size_bp;
    let checkpoint = chromosome_index
        .bins
        .iter()
        .filter(|row| row.bin_index >= first_bin && row.bin_index <= last_bin)
        .min_by_key(|row| row.scan_start_byte);
    let Some(checkpoint) = checkpoint else {
        return Ok(vec![]);
    };
    let accepted_types = requested_types
        .iter()
        .map(|value| value.trim().to_ascii_lowercase())
        .filter(|value| !value.is_empty())
        .collect::<BTreeSet<_>>();
    let mut file = File::open(intervals_path).map_err(|error| {
        format!(
            "Could not open Ensembl Regulation intervals '{}': {error}",
            intervals_path.display()
        )
    })?;
    file.seek(SeekFrom::Start(checkpoint.scan_start_byte))
        .map_err(|error| {
            format!(
                "Could not seek Ensembl Regulation intervals '{}': {error}",
                intervals_path.display()
            )
        })?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let mut line_number = checkpoint.scan_start_line_number;
    let mut rows = Vec::new();
    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line).map_err(|error| {
            format!(
                "Could not query Ensembl Regulation intervals '{}': {error}",
                intervals_path.display()
            )
        })?;
        if bytes_read == 0 {
            break;
        }
        let Some(interval) = parse_interval_line(&line, line_number)? else {
            line_number = line_number.saturating_add(1);
            continue;
        };
        if interval.chromosome != canonical {
            break;
        }
        if interval.start_0based >= end_0based_exclusive {
            break;
        }
        if interval.end_0based_exclusive > start_0based
            && (accepted_types.is_empty()
                || accepted_types.contains(&interval.feature_type.to_ascii_lowercase()))
        {
            rows.push(interval);
        }
        line_number = line_number.saturating_add(1);
    }
    rows.sort_by(|left, right| {
        left.start_0based
            .cmp(&right.start_0based)
            .then(left.end_0based_exclusive.cmp(&right.end_0based_exclusive))
            .then(left.feature_id.cmp(&right.feature_id))
    });
    Ok(rows)
}

pub fn genome_id_matches_source(
    genome_id: &str,
    source: &EnsemblRegulationSourceDescriptor,
) -> bool {
    let tokens = genome_id
        .split(|ch: char| !ch.is_ascii_alphanumeric())
        .filter(|value| !value.is_empty())
        .collect::<Vec<_>>();
    source.assembly_aliases.iter().any(|alias| {
        alias.eq_ignore_ascii_case(genome_id)
            || tokens.iter().any(|token| token.eq_ignore_ascii_case(alias))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn fixture_response() -> ApiResponse {
        serde_json::from_str(
            r#"{
              "release":"2026-08",
              "regions":[{
                "region":{"name":"1","length":100000,"coordinate_system":"chromosome"},
                "features":[
                  {"id":"ENSR_TEST_PROMOTER","feature_type":"promoter","start":101,"end":180,"strand":"independent","extended_start":91,"extended_end":190,"associated_genes":["ENSG_TEST"],"associated_gene_names":["GENE1"]},
                  {"id":"ENSR_TEST_CTCF","feature_type":"ctcf","start":65521,"end":65580,"strand":"forward","associated_genes":[],"associated_gene_names":[]}
                ]
              }],
              "limit":2,
              "total_count":2
            }"#,
        )
        .expect("fixture response")
    }

    fn write_fixture(path: &Path) {
        let source =
            source_descriptor(ENSEMBL_REGULATION_2026_08_GRCH38_SOURCE_ID).expect("human source");
        let rows = intervals_from_api_response(&fixture_response(), &source).expect("rows");
        let file = File::create(path).expect("interval fixture");
        let mut writer = BufWriter::new(file);
        write_intervals_header(&mut writer, &source).expect("header");
        for row in &rows {
            write_interval(&mut writer, row).expect("row");
        }
        writer.flush().expect("flush");
    }

    #[test]
    fn catalog_keeps_ensembl_assemblies_distinct_from_screen_mouse() {
        let human =
            source_descriptor(ENSEMBL_REGULATION_2026_08_GRCH38_SOURCE_ID).expect("human source");
        let mouse =
            source_descriptor(ENSEMBL_REGULATION_2026_08_GRCM39_SOURCE_ID).expect("mouse source");
        assert_eq!(human.gene_annotation_release, "50");
        assert_eq!(mouse.assembly_name, "GRCm39");
        assert!(genome_id_matches_source("mm39_ensembl", &mouse));
        assert!(!genome_id_matches_source("mm10", &mouse));
        assert!(!genome_id_matches_source("hg38", &mouse));
        assert_eq!(human.browser_species_slug, "homo_sapiens");
        assert_eq!(mouse.browser_species_slug, "mus_musculus");
    }

    #[test]
    fn feature_urls_use_only_the_pinned_provider_contract() {
        let human =
            source_descriptor(ENSEMBL_REGULATION_2026_08_GRCH38_SOURCE_ID).expect("human source");
        assert_eq!(
            feature_page_url(&human, "ENSR1_958").expect("feature URL"),
            "https://regulation.ensembl.org/regulatory_features/homo_sapiens/ENSR1_958"
        );
        assert!(feature_page_url(&human, "https://example.test/x").is_err());
        let mut tampered = human;
        tampered.feature_page_url_template = "https://example.test/{FEATURE_ID}".to_string();
        assert!(feature_page_url(&tampered, "ENSR1_958").is_err());
    }

    #[test]
    fn api_coordinates_and_provider_metadata_are_preserved() {
        let source =
            source_descriptor(ENSEMBL_REGULATION_2026_08_GRCH38_SOURCE_ID).expect("human source");
        let rows = intervals_from_api_response(&fixture_response(), &source).expect("rows");
        assert_eq!(rows[0].start_0based, 100);
        assert_eq!(rows[0].end_0based_exclusive, 180);
        assert_eq!(rows[0].extended_start_0based, Some(90));
        assert_eq!(rows[0].extended_end_0based_exclusive, Some(190));
        assert_eq!(rows[0].associated_gene_names, ["GENE1"]);
        assert_eq!(rows[1].strand, Some('+'));
    }

    #[test]
    fn sparse_index_queries_aliases_and_verifies_content() {
        let root = tempdir().expect("tempdir");
        let intervals = root.path().join("regulation.tsv");
        let index_path = root.path().join("regulation.index.json");
        write_fixture(&intervals);
        let index = prepare_interval_index(
            intervals.to_string_lossy().as_ref(),
            index_path.to_string_lossy().as_ref(),
            source_descriptor(ENSEMBL_REGULATION_2026_08_GRCH38_SOURCE_ID).expect("source"),
        )
        .expect("index");
        assert_eq!(index.row_count, 2);
        validate_intervals_identity(&intervals, &index).expect("identity");
        let rows = query_overlaps(&index, &intervals, "chr1", 65_570, 65_575, &[]).expect("query");
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].feature_id, "ENSR_TEST_CTCF");
    }

    #[test]
    fn legacy_index_hydrates_empty_browser_contract_from_pinned_source() {
        let root = tempdir().expect("tempdir");
        let intervals = root.path().join("regulation.tsv");
        let index_path = root.path().join("regulation.index.json");
        write_fixture(&intervals);
        let index = prepare_interval_index(
            intervals.to_string_lossy().as_ref(),
            index_path.to_string_lossy().as_ref(),
            source_descriptor(ENSEMBL_REGULATION_2026_08_GRCH38_SOURCE_ID).expect("source"),
        )
        .expect("index");
        let mut value = serde_json::to_value(index).expect("serialize index");
        let source = value["source"].as_object_mut().expect("source object");
        source.remove("browser_species_slug");
        source.remove("feature_page_url_template");
        fs::write(
            &index_path,
            serde_json::to_vec_pretty(&value).expect("serialize legacy index"),
        )
        .expect("write legacy index");
        let hydrated =
            read_interval_index(index_path.to_string_lossy().as_ref()).expect("read legacy index");
        assert_eq!(hydrated.source.browser_species_slug, "homo_sapiens");
        assert!(
            hydrated
                .source
                .feature_page_url_template
                .contains("{FEATURE_ID}")
        );
    }

    #[test]
    fn same_size_interval_content_replacement_fails_digest_validation() {
        let root = tempdir().expect("tempdir");
        let intervals = root.path().join("regulation.tsv");
        let index_path = root.path().join("regulation.index.json");
        write_fixture(&intervals);
        let index = prepare_interval_index(
            intervals.to_string_lossy().as_ref(),
            index_path.to_string_lossy().as_ref(),
            source_descriptor(ENSEMBL_REGULATION_2026_08_GRCH38_SOURCE_ID).expect("source"),
        )
        .expect("index");
        let mut bytes = fs::read(&intervals).expect("read intervals");
        let marker = b"ENSR_TEST_PROMOTER";
        let marker_start = bytes
            .windows(marker.len())
            .position(|window| window == marker)
            .expect("fixture feature id");
        bytes[marker_start + marker.len() - 1] = b'X';
        fs::write(&intervals, bytes).expect("replace interval content in place");

        let error = validate_intervals_identity(&intervals, &index)
            .expect_err("same-size content replacement must fail its digest binding");
        assert!(error.contains("content digest"));
    }

    #[test]
    fn release_mismatch_fails_closed() {
        let source =
            source_descriptor(ENSEMBL_REGULATION_2026_08_GRCH38_SOURCE_ID).expect("human source");
        let mut response = fixture_response();
        response.release = "2025-05".to_string();
        let error = intervals_from_api_response(&response, &source).expect_err("release mismatch");
        assert!(error.contains("stale or redirected"));
    }
}
