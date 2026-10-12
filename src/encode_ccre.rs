//! Optional ENCODE SCREEN Registry cCRE source, index, and overlap helpers.
//!
//! The large BED payload stays external to project state. A compact sparse
//! byte-offset index is derived from the sorted BED file and bound to its exact
//! SHA-256 digest, allowing explicit local interval queries without scanning or
//! loading the resource during normal GENtle startup.

use crate::digest_utils::sha256_file_hex;
use gentle_protocol::{
    ENCODE_CCRE_INTERVAL_INDEX_SCHEMA, ENCODE_CCRE_SOURCE_CATALOG_SCHEMA, EncodeCcreBinCheckpoint,
    EncodeCcreChromosomeIndex, EncodeCcreInterval, EncodeCcreIntervalIndex,
    EncodeCcreSourceCatalog, EncodeCcreSourceDescriptor,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
};

pub const DEFAULT_ENCODE_CCRE_BIN_SIZE_BP: u64 = 65_536;
pub const SCREEN_V4_GRCH38_ELS_SOURCE_ID: &str = "screen_registry_v4_grch38_els";
pub const SCREEN_V4_MM10_ELS_SOURCE_ID: &str = "screen_registry_v4_mm10_els";
pub const SCREEN_DOWNLOAD_PAGE_URL: &str = "https://screen.wenglab.org/downloads";
pub const SCREEN_V4_PUBLICATION_URL: &str = "https://www.nature.com/articles/s41586-025-09909-9";

fn els_source(
    source_id: &str,
    species: &str,
    taxon_id: u32,
    assembly_name: &str,
    assembly_aliases: &[&str],
    dhs_accession_prefix: &str,
    ccre_accession_prefix: &str,
    source_url: &str,
) -> EncodeCcreSourceDescriptor {
    EncodeCcreSourceDescriptor {
        source_id: source_id.to_string(),
        provider: "ENCODE SCREEN".to_string(),
        registry_version: "4".to_string(),
        species_scientific_name: species.to_string(),
        taxon_id,
        assembly_name: assembly_name.to_string(),
        assembly_aliases: assembly_aliases
            .iter()
            .map(|value| (*value).to_string())
            .collect(),
        subset_id: "ELS".to_string(),
        subset_label: "All candidate enhancer-like signatures (pELS and dELS)".to_string(),
        primary_classes: vec!["pELS".to_string(), "dELS".to_string()],
        dhs_accession_prefix: dhs_accession_prefix.to_string(),
        ccre_accession_prefix: ccre_accession_prefix.to_string(),
        source_url: source_url.to_string(),
        download_page_url: SCREEN_DOWNLOAD_PAGE_URL.to_string(),
        publication_url: SCREEN_V4_PUBLICATION_URL.to_string(),
        coordinate_system: "BED 0-based half-open".to_string(),
        field_order: vec![
            "chromosome".to_string(),
            "start_0based".to_string(),
            "end_0based_exclusive".to_string(),
            "dhs_accession".to_string(),
            "ccre_accession".to_string(),
            "ccre_class".to_string(),
        ],
        scope_note: "Enhancer-like cCRE subset only; this is not the complete Registry V4 cCRE catalog. cCRE class is biochemical candidate evidence, not proof of enhancer activity or target-gene regulation."
            .to_string(),
    }
}

pub fn source_catalog() -> EncodeCcreSourceCatalog {
    EncodeCcreSourceCatalog {
        schema: ENCODE_CCRE_SOURCE_CATALOG_SCHEMA.to_string(),
        registry_version: "4".to_string(),
        sources: vec![
            els_source(
                SCREEN_V4_GRCH38_ELS_SOURCE_ID,
                "Homo sapiens",
                9606,
                "GRCh38",
                &["GRCh38", "hg38"],
                "EH38D",
                "EH38E",
                "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.ELS.bed",
            ),
            els_source(
                SCREEN_V4_MM10_ELS_SOURCE_ID,
                "Mus musculus",
                10090,
                "GRCm38",
                &["GRCm38", "mm10"],
                "EM10D",
                "EM10E",
                "https://downloads.wenglab.org/Registry-V4/mm10-cCREs.ELS.bed",
            ),
        ],
        notes: vec![
            "Registry V4 is the default; older Registry V3 URLs are not silently substituted."
                .to_string(),
            "The human and mouse ELS files are separate species- and assembly-bound resources."
                .to_string(),
            "No resource is downloaded, opened, or indexed until an explicit resource command is run."
                .to_string(),
        ],
    }
}

pub fn source_descriptor(source_id: &str) -> Result<EncodeCcreSourceDescriptor, String> {
    let source_id = source_id.trim();
    source_catalog()
        .sources
        .into_iter()
        .find(|row| row.source_id == source_id)
        .ok_or_else(|| {
            format!(
                "Unknown ENCODE SCREEN cCRE source '{source_id}'; available source ids: {}",
                source_catalog()
                    .sources
                    .iter()
                    .map(|row| row.source_id.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        })
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

pub fn default_bed_path(source: &EncodeCcreSourceDescriptor) -> String {
    format!(
        "data/resources/screen.registry-v{}.{}.{}.bed",
        path_token(&source.registry_version),
        path_token(
            source
                .assembly_aliases
                .last()
                .unwrap_or(&source.assembly_name)
        ),
        path_token(&source.subset_id)
    )
}

pub fn default_index_path(source: &EncodeCcreSourceDescriptor) -> String {
    format!("{}.interval-index.json", default_bed_path(source))
}

pub fn parse_bed_line(
    line: &str,
    source_line_number: u64,
) -> Result<Option<EncodeCcreInterval>, String> {
    let trimmed = line.trim_end_matches(['\r', '\n']);
    if trimmed.trim().is_empty()
        || trimmed.starts_with('#')
        || trimmed.starts_with("track ")
        || trimmed.starts_with("browser ")
    {
        return Ok(None);
    }
    let fields = trimmed.split('\t').collect::<Vec<_>>();
    if fields.len() != 6 {
        return Err(format!(
            "SCREEN cCRE BED line {source_line_number} has {} fields; expected exactly 6",
            fields.len()
        ));
    }
    let chromosome = fields[0].trim();
    let start_0based = fields[1].parse::<u64>().map_err(|error| {
        format!("Invalid SCREEN cCRE start on line {source_line_number}: {error}")
    })?;
    let end_0based_exclusive = fields[2].parse::<u64>().map_err(|error| {
        format!("Invalid SCREEN cCRE end on line {source_line_number}: {error}")
    })?;
    if chromosome.is_empty() || end_0based_exclusive <= start_0based {
        return Err(format!(
            "Invalid SCREEN cCRE interval on line {source_line_number}: '{chromosome}:{}..{}'",
            start_0based, end_0based_exclusive
        ));
    }
    for (label, value) in [
        ("DHS accession", fields[3]),
        ("cCRE accession", fields[4]),
        ("cCRE class", fields[5]),
    ] {
        if value.trim().is_empty() {
            return Err(format!(
                "SCREEN cCRE {label} is blank on line {source_line_number}"
            ));
        }
    }
    Ok(Some(EncodeCcreInterval {
        chromosome: chromosome.to_string(),
        start_0based,
        end_0based_exclusive,
        dhs_accession: fields[3].trim().to_string(),
        ccre_accession: fields[4].trim().to_string(),
        ccre_class: fields[5].trim().to_string(),
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

fn write_index(path: &str, index: &EncodeCcreIntervalIndex) -> Result<(), String> {
    let output = Path::new(path);
    if let Some(parent) = output
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent).map_err(|error| {
            format!(
                "Could not create ENCODE cCRE index directory '{}': {error}",
                parent.display()
            )
        })?;
    }
    let file = File::create(output)
        .map_err(|error| format!("Could not create ENCODE cCRE index '{path}': {error}"))?;
    let mut writer = BufWriter::new(file);
    serde_json::to_writer_pretty(&mut writer, index)
        .map_err(|error| format!("Could not serialize ENCODE cCRE index '{path}': {error}"))?;
    writer
        .write_all(b"\n")
        .and_then(|_| writer.flush())
        .map_err(|error| format!("Could not finish ENCODE cCRE index '{path}': {error}"))
}

pub fn prepare_interval_index(
    bed_path: &str,
    output_path: &str,
    source: EncodeCcreSourceDescriptor,
) -> Result<EncodeCcreIntervalIndex, String> {
    if source.registry_version != "4" {
        return Err(format!(
            "ENCODE cCRE index preparation currently requires Registry V4, got '{}'",
            source.registry_version
        ));
    }
    let file = File::open(bed_path)
        .map_err(|error| format!("Could not open ENCODE cCRE BED '{bed_path}': {error}"))?;
    let bed_size_bytes = file
        .metadata()
        .map_err(|error| format!("Could not inspect ENCODE cCRE BED '{bed_path}': {error}"))?
        .len();
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let mut byte_offset = 0u64;
    let mut source_line_number = 0u64;
    let mut row_count = 0u64;
    let mut class_counts = BTreeMap::<String, u64>::new();
    let mut chromosome_bins = BTreeMap::<String, BTreeMap<u64, (u64, u64)>>::new();
    let mut chromosome_meta = BTreeMap::<String, (u64, u64, u64)>::new();
    let mut seen_chromosomes = BTreeSet::<String>::new();
    let mut current_chromosome: Option<String> = None;
    let mut previous_start = 0u64;
    let accepted_classes = source
        .primary_classes
        .iter()
        .map(|value| value.to_ascii_lowercase())
        .collect::<BTreeSet<_>>();

    loop {
        line.clear();
        let bytes_read = reader
            .read_line(&mut line)
            .map_err(|error| format!("Could not read ENCODE cCRE BED '{bed_path}': {error}"))?;
        if bytes_read == 0 {
            break;
        }
        source_line_number = source_line_number.saturating_add(1);
        let next_byte_offset = byte_offset.saturating_add(bytes_read as u64);
        let Some(interval) = parse_bed_line(&line, source_line_number)? else {
            byte_offset = next_byte_offset;
            continue;
        };
        if !accepted_classes.is_empty()
            && !accepted_classes.contains(&interval.ccre_class.to_ascii_lowercase())
        {
            return Err(format!(
                "SCREEN cCRE BED line {} has class '{}' outside source '{}' declared classes {}",
                source_line_number,
                interval.ccre_class,
                source.source_id,
                source.primary_classes.join(", ")
            ));
        }
        if !source.dhs_accession_prefix.is_empty()
            && !interval
                .dhs_accession
                .starts_with(&source.dhs_accession_prefix)
        {
            return Err(format!(
                "SCREEN cCRE BED line {} has DHS accession '{}' outside source '{}' expected prefix '{}'; verify species and assembly",
                source_line_number,
                interval.dhs_accession,
                source.source_id,
                source.dhs_accession_prefix
            ));
        }
        if !source.ccre_accession_prefix.is_empty()
            && !interval
                .ccre_accession
                .starts_with(&source.ccre_accession_prefix)
        {
            return Err(format!(
                "SCREEN cCRE BED line {} has cCRE accession '{}' outside source '{}' expected prefix '{}'; verify species and assembly",
                source_line_number,
                interval.ccre_accession,
                source.source_id,
                source.ccre_accession_prefix
            ));
        }
        match current_chromosome.as_deref() {
            Some(current) if current == interval.chromosome => {
                if interval.start_0based < previous_start {
                    return Err(format!(
                        "SCREEN cCRE BED is not start-sorted on chromosome '{}' at line {}",
                        interval.chromosome, source_line_number
                    ));
                }
            }
            Some(_) => {
                if seen_chromosomes.contains(&interval.chromosome) {
                    return Err(format!(
                        "SCREEN cCRE BED chromosome '{}' is split across non-contiguous blocks (line {})",
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
        let first_bin = interval.start_0based / DEFAULT_ENCODE_CCRE_BIN_SIZE_BP;
        let last_bin =
            interval.end_0based_exclusive.saturating_sub(1) / DEFAULT_ENCODE_CCRE_BIN_SIZE_BP;
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
        *class_counts.entry(interval.ccre_class).or_insert(0) += 1;
        row_count = row_count.saturating_add(1);
        byte_offset = next_byte_offset;
    }
    if row_count == 0 {
        return Err(format!(
            "ENCODE cCRE BED '{bed_path}' contains no data rows"
        ));
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
            EncodeCcreChromosomeIndex {
                chromosome,
                row_count: chromosome_row_count,
                first_row_byte,
                end_byte_exclusive,
                bins: bins
                    .into_iter()
                    .map(|(bin_index, (scan_start_byte, scan_start_line_number))| {
                        EncodeCcreBinCheckpoint {
                            bin_index,
                            scan_start_byte,
                            scan_start_line_number,
                        }
                    })
                    .collect(),
            },
        );
    }
    let bed_sha256 = format!(
        "sha256:{}",
        sha256_file_hex(Path::new(bed_path))
            .map_err(|error| { format!("Could not hash ENCODE cCRE BED '{bed_path}': {error}") })?
    );
    let index = EncodeCcreIntervalIndex {
        schema: ENCODE_CCRE_INTERVAL_INDEX_SCHEMA.to_string(),
        source,
        bed_path_hint: bed_path.to_string(),
        bed_file_name: Path::new(bed_path)
            .file_name()
            .and_then(|value| value.to_str())
            .unwrap_or_default()
            .to_string(),
        bed_sha256,
        bed_size_bytes,
        row_count,
        chromosome_count: chromosomes.len(),
        class_counts,
        bin_size_bp: DEFAULT_ENCODE_CCRE_BIN_SIZE_BP,
        chromosomes,
        chromosome_aliases: chromosome_alias_map,
        warnings: vec![],
    };
    write_index(output_path, &index)?;
    Ok(index)
}

pub fn read_interval_index(path: &str) -> Result<EncodeCcreIntervalIndex, String> {
    let text = fs::read_to_string(path)
        .map_err(|error| format!("Could not read ENCODE cCRE index '{path}': {error}"))?;
    let index = serde_json::from_str::<EncodeCcreIntervalIndex>(&text)
        .map_err(|error| format!("Invalid ENCODE cCRE index '{path}': {error}"))?;
    if index.schema != ENCODE_CCRE_INTERVAL_INDEX_SCHEMA {
        return Err(format!(
            "Unsupported ENCODE cCRE index schema '{}' in '{}'; expected '{}'",
            index.schema, path, ENCODE_CCRE_INTERVAL_INDEX_SCHEMA
        ));
    }
    if index.bin_size_bp == 0 || index.bed_sha256.trim().is_empty() {
        return Err(format!(
            "ENCODE cCRE index '{path}' lacks a valid bin size or BED content digest"
        ));
    }
    let catalog_source = source_descriptor(&index.source.source_id).map_err(|error| {
        format!("ENCODE cCRE index '{path}' does not identify a supported source: {error}")
    })?;
    if index.source != catalog_source {
        return Err(format!(
            "ENCODE cCRE index '{path}' source descriptor does not match GENtle's '{}' Registry V4 catalog entry; rebuild the index before use",
            index.source.source_id
        ));
    }
    Ok(index)
}

pub fn resolve_bed_path(
    index_path: &str,
    index: &EncodeCcreIntervalIndex,
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
            "Explicit ENCODE cCRE BED override '{}' does not exist",
            path.display()
        ));
    }
    if !index.bed_file_name.is_empty() {
        let sibling = Path::new(index_path)
            .parent()
            .unwrap_or_else(|| Path::new(""))
            .join(&index.bed_file_name);
        if sibling.is_file() {
            return Ok(sibling);
        }
    }
    let hint = PathBuf::from(&index.bed_path_hint);
    if hint.is_file() {
        return Ok(hint);
    }
    Err(format!(
        "Could not locate ENCODE cCRE BED for index '{}'; expected sibling '{}' or recorded path '{}'. Pass --bed PATH after relocation.",
        index_path, index.bed_file_name, index.bed_path_hint
    ))
}

fn resolve_chromosome<'a>(index: &'a EncodeCcreIntervalIndex, requested: &str) -> Option<&'a str> {
    if let Some(canonical) = index.chromosome_aliases.get(requested) {
        return Some(canonical.as_str());
    }
    index
        .chromosome_aliases
        .iter()
        .find(|(alias, _)| alias.eq_ignore_ascii_case(requested))
        .map(|(_, canonical)| canonical.as_str())
}

pub fn validate_bed_identity(path: &Path, index: &EncodeCcreIntervalIndex) -> Result<(), String> {
    let metadata = fs::metadata(path).map_err(|error| {
        format!(
            "Could not inspect ENCODE cCRE BED '{}': {error}",
            path.display()
        )
    })?;
    if metadata.len() != index.bed_size_bytes {
        return Err(format!(
            "ENCODE cCRE BED '{}' has {} bytes, but the index binds {} bytes; rebuild the index",
            path.display(),
            metadata.len(),
            index.bed_size_bytes
        ));
    }
    let observed = format!(
        "sha256:{}",
        sha256_file_hex(path).map_err(|error| {
            format!(
                "Could not hash ENCODE cCRE BED '{}': {error}",
                path.display()
            )
        })?
    );
    if observed != index.bed_sha256 {
        return Err(format!(
            "ENCODE cCRE BED '{}' content digest '{}' does not match index digest '{}'; rebuild the index",
            path.display(),
            observed,
            index.bed_sha256
        ));
    }
    Ok(())
}

pub fn query_overlaps(
    index: &EncodeCcreIntervalIndex,
    bed_path: &Path,
    chromosome: &str,
    start_0based: u64,
    end_0based_exclusive: u64,
    requested_classes: &[String],
) -> Result<Vec<EncodeCcreInterval>, String> {
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
    let accepted_classes = requested_classes
        .iter()
        .map(|value| value.trim().to_ascii_lowercase())
        .filter(|value| !value.is_empty())
        .collect::<BTreeSet<_>>();
    let mut file = File::open(bed_path).map_err(|error| {
        format!(
            "Could not open ENCODE cCRE BED '{}': {error}",
            bed_path.display()
        )
    })?;
    file.seek(SeekFrom::Start(checkpoint.scan_start_byte))
        .map_err(|error| {
            format!(
                "Could not seek ENCODE cCRE BED '{}': {error}",
                bed_path.display()
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
                "Could not query ENCODE cCRE BED '{}': {error}",
                bed_path.display()
            )
        })?;
        if bytes_read == 0 {
            break;
        }
        let Some(interval) = parse_bed_line(&line, line_number)? else {
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
            && (accepted_classes.is_empty()
                || accepted_classes.contains(&interval.ccre_class.to_ascii_lowercase()))
        {
            rows.push(interval);
        }
        line_number = line_number.saturating_add(1);
    }
    rows.sort_by(|left, right| {
        left.start_0based
            .cmp(&right.start_0based)
            .then(left.end_0based_exclusive.cmp(&right.end_0based_exclusive))
            .then(left.ccre_accession.cmp(&right.ccre_accession))
            .then(left.dhs_accession.cmp(&right.dhs_accession))
    });
    Ok(rows)
}

pub fn genome_id_matches_source(genome_id: &str, source: &EncodeCcreSourceDescriptor) -> bool {
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

    const FIXTURE: &str = concat!(
        "chr1\t100\t180\tEH38D0000001\tEH38E0000001\tpELS\n",
        "chr1\t65520\t65580\tEH38D0000002\tEH38E0000002\tdELS\n",
        "chr1\t70000\t70100\tEH38D0000003\tEH38E0000003\tdELS\n",
        "chr2\t10\t20\tEH38D0000004\tEH38E0000004\tpELS\n",
    );

    #[test]
    fn catalog_exposes_separate_human_and_mouse_v4_els_sources() {
        let catalog = source_catalog();
        assert_eq!(catalog.registry_version, "4");
        assert_eq!(catalog.sources.len(), 2);
        let human = source_descriptor(SCREEN_V4_GRCH38_ELS_SOURCE_ID).expect("human source");
        assert_eq!(human.taxon_id, 9606);
        assert_eq!(human.assembly_name, "GRCh38");
        assert_eq!(human.primary_classes, ["pELS", "dELS"]);
        assert_eq!(human.dhs_accession_prefix, "EH38D");
        assert_eq!(human.ccre_accession_prefix, "EH38E");
        assert_eq!(
            human.source_url,
            "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.ELS.bed"
        );
        let mouse = source_descriptor(SCREEN_V4_MM10_ELS_SOURCE_ID).expect("mouse source");
        assert_eq!(mouse.taxon_id, 10090);
        assert_eq!(mouse.assembly_name, "GRCm38");
        assert_eq!(mouse.dhs_accession_prefix, "EM10D");
        assert_eq!(mouse.ccre_accession_prefix, "EM10E");
        assert_eq!(
            mouse.source_url,
            "https://downloads.wenglab.org/Registry-V4/mm10-cCREs.ELS.bed"
        );
    }

    #[test]
    fn sparse_index_finds_an_interval_crossing_a_bin_boundary() {
        let root = tempdir().expect("tempdir");
        let bed = root.path().join("ccres.bed");
        let index_path = root.path().join("ccres.index.json");
        fs::write(&bed, FIXTURE).expect("fixture BED");
        let source = source_descriptor(SCREEN_V4_GRCH38_ELS_SOURCE_ID).expect("source");
        let index = prepare_interval_index(
            bed.to_string_lossy().as_ref(),
            index_path.to_string_lossy().as_ref(),
            source,
        )
        .expect("prepare index");
        assert_eq!(index.row_count, 4);
        assert_eq!(index.class_counts["pELS"], 2);
        assert_eq!(index.class_counts["dELS"], 2);
        validate_bed_identity(&bed, &index).expect("identity");
        let rows = query_overlaps(&index, &bed, "1", 65_570, 65_575, &[]).expect("boundary query");
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].ccre_accession, "EH38E0000002");
    }

    #[test]
    fn query_class_filter_and_replaced_payload_fail_closed() {
        let root = tempdir().expect("tempdir");
        let bed = root.path().join("ccres.bed");
        let index_path = root.path().join("ccres.index.json");
        fs::write(&bed, FIXTURE).expect("fixture BED");
        let index = prepare_interval_index(
            bed.to_string_lossy().as_ref(),
            index_path.to_string_lossy().as_ref(),
            source_descriptor(SCREEN_V4_GRCH38_ELS_SOURCE_ID).expect("source"),
        )
        .expect("prepare index");
        let rows = query_overlaps(&index, &bed, "chr1", 0, 80_000, &["pELS".to_string()])
            .expect("filtered query");
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].ccre_class, "pELS");
        let mut replacement = FIXTURE.as_bytes().to_vec();
        replacement[5] = b'2';
        fs::write(&bed, replacement).expect("replace BED");
        let error = validate_bed_identity(&bed, &index).expect_err("digest mismatch");
        assert!(error.contains("content digest"));
    }

    #[test]
    fn assembly_matching_keeps_human_and_mouse_separate() {
        let human = source_descriptor(SCREEN_V4_GRCH38_ELS_SOURCE_ID).expect("human");
        let mouse = source_descriptor(SCREEN_V4_MM10_ELS_SOURCE_ID).expect("mouse");
        assert!(genome_id_matches_source("grch38_ensembl116", &human));
        assert!(!genome_id_matches_source("mm10_ensembl116", &human));
        assert!(genome_id_matches_source("mm10_ensembl116", &mouse));
        assert!(!genome_id_matches_source("hg38", &mouse));
    }

    #[test]
    fn index_rejects_a_bed_from_the_wrong_species() {
        let root = tempdir().expect("tempdir");
        let bed = root.path().join("human-labelled-as-mouse.bed");
        let index_path = root.path().join("ccres.index.json");
        fs::write(&bed, FIXTURE).expect("fixture BED");
        let error = prepare_interval_index(
            bed.to_string_lossy().as_ref(),
            index_path.to_string_lossy().as_ref(),
            source_descriptor(SCREEN_V4_MM10_ELS_SOURCE_ID).expect("mouse source"),
        )
        .expect_err("human accessions must not satisfy mouse source identity");
        assert!(error.contains("expected prefix 'EM10D'"));
        assert!(error.contains("verify species and assembly"));
    }

    #[test]
    fn index_load_rejects_a_tampered_species_descriptor() {
        let root = tempdir().expect("tempdir");
        let bed = root.path().join("ccres.bed");
        let index_path = root.path().join("ccres.index.json");
        fs::write(&bed, FIXTURE).expect("fixture BED");
        let mut index = prepare_interval_index(
            bed.to_string_lossy().as_ref(),
            index_path.to_string_lossy().as_ref(),
            source_descriptor(SCREEN_V4_GRCH38_ELS_SOURCE_ID).expect("human source"),
        )
        .expect("prepare index");
        index.source.species_scientific_name = "Mus musculus".to_string();
        write_index(index_path.to_string_lossy().as_ref(), &index).expect("tampered index");

        let error = read_interval_index(index_path.to_string_lossy().as_ref())
            .expect_err("tampered species descriptor must fail closed");
        assert!(error.contains("source descriptor does not match"));
        assert!(error.contains(SCREEN_V4_GRCH38_ELS_SOURCE_ID));
    }
}
