//! microRNA seed-site scanning over annotated DNA sequence records.
//!
//! The scanner reports sequence evidence only: exact canonical seed matches
//! are candidates, not functional validation.

use crate::{
    dna_sequence::DNAsequence,
    feature_location::{collect_location_ranges_usize, feature_is_reverse},
};
use gb_io::seq::Feature;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const MIRNA_TARGET_SCAN_SCHEMA: &str = "gentle.mirna_target_scan.v1";
pub const MIRNA_SEED_EXPLANATION_SCHEMA: &str = "gentle.mirna_seed_explanation.v1";
pub const MIRNA_CATALOG_RECORD_SCHEMA: &str = "gentle.mirna_catalog_record.v1";
const DEFAULT_CONTEXT_BP: usize = 20;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum MirnaSeedClass {
    #[serde(rename = "8mer")]
    EightMer,
    #[serde(rename = "7mer-m8")]
    SevenMerM8,
    #[serde(rename = "7mer-A1")]
    SevenMerA1,
    #[serde(rename = "6mer")]
    SixMer,
}

impl MirnaSeedClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::EightMer => "8mer",
            Self::SevenMerM8 => "7mer-m8",
            Self::SevenMerA1 => "7mer-A1",
            Self::SixMer => "6mer",
        }
    }

    pub fn parse(raw: &str) -> Result<Self, String> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "8mer" | "8" => Ok(Self::EightMer),
            "7mer-m8" | "7merm8" | "m8" => Ok(Self::SevenMerM8),
            "7mer-a1" | "7mera1" | "a1" => Ok(Self::SevenMerA1),
            "6mer" | "6" => Ok(Self::SixMer),
            other => Err(format!("Unknown microRNA seed class '{other}'")),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum MirnaRegionClass {
    #[serde(rename = "3_prime_utr")]
    ThreePrimeUtr,
    #[serde(rename = "coding_exon")]
    CodingExon,
    #[serde(rename = "noncoding_exon")]
    NoncodingExon,
    #[serde(rename = "intron")]
    Intron,
    #[serde(rename = "exon_intron_boundary")]
    ExonIntronBoundary,
    #[serde(rename = "intron_exon_boundary")]
    IntronExonBoundary,
    #[serde(rename = "whole_transcript")]
    WholeTranscript,
    #[serde(rename = "whole_gene")]
    WholeGene,
}

impl MirnaRegionClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ThreePrimeUtr => "3_prime_utr",
            Self::CodingExon => "coding_exon",
            Self::NoncodingExon => "noncoding_exon",
            Self::Intron => "intron",
            Self::ExonIntronBoundary => "exon_intron_boundary",
            Self::IntronExonBoundary => "intron_exon_boundary",
            Self::WholeTranscript => "whole_transcript",
            Self::WholeGene => "whole_gene",
        }
    }

    pub fn parse(raw: &str) -> Result<Vec<Self>, String> {
        let value = raw.trim().to_ascii_lowercase();
        match value.as_str() {
            "3utr" | "3'utr" | "3_prime_utr" | "three_prime_utr" => Ok(vec![Self::ThreePrimeUtr]),
            "coding_exon" | "coding-exon" => Ok(vec![Self::CodingExon]),
            "noncoding_exon" | "noncoding-exon" | "utr_exon" | "utr-exon" => {
                Ok(vec![Self::NoncodingExon])
            }
            "exon" | "exons" => Ok(vec![Self::CodingExon, Self::NoncodingExon]),
            "intron" | "introns" => Ok(vec![Self::Intron]),
            "exon_intron_boundary" | "exon-intron-boundary" => Ok(vec![Self::ExonIntronBoundary]),
            "intron_exon_boundary" | "intron-exon-boundary" => Ok(vec![Self::IntronExonBoundary]),
            "boundary" | "boundaries" => {
                Ok(vec![Self::ExonIntronBoundary, Self::IntronExonBoundary])
            }
            "whole_transcript" | "whole-transcript" => Ok(vec![Self::WholeTranscript]),
            "whole_gene" | "whole-gene" => Ok(vec![Self::WholeGene]),
            other => Err(format!("Unknown microRNA scan region '{other}'")),
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MirnaCatalogRecord {
    pub schema: String,
    pub id: String,
    pub aliases: Vec<String>,
    pub mature_sequence: String,
    pub accession: Option<String>,
    pub source: String,
    pub ncbi_gene_id: Option<String>,
    pub notes: Vec<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MirnaQueryRecord {
    pub id: String,
    pub mature_sequence: String,
    pub accession: Option<String>,
    pub source: String,
    pub aliases: Vec<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MirnaSeedMotif {
    pub seed_class: MirnaSeedClass,
    pub mirna_positions_1based: String,
    pub mirna_seed_5p: String,
    pub target_motif: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MirnaTargetScanRequest {
    pub mirna: String,
    pub mature_sequence: Option<String>,
    pub regions: Vec<MirnaRegionClass>,
    pub seed_classes: Vec<MirnaSeedClass>,
    pub boundary_flank_bp: usize,
    pub transcript_filter: Option<String>,
    pub species_note: Option<String>,
    pub evidence_notes: Vec<String>,
}

impl MirnaTargetScanRequest {
    pub fn with_defaults(mirna: impl Into<String>) -> Self {
        Self {
            mirna: mirna.into(),
            mature_sequence: None,
            regions: default_scan_regions(),
            seed_classes: default_seed_classes(),
            boundary_flank_bp: 25,
            transcript_filter: None,
            species_note: None,
            evidence_notes: vec![],
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MirnaTargetHit {
    pub region_class: MirnaRegionClass,
    pub transcript_id: Option<String>,
    pub feature_id: Option<String>,
    pub parent_feature: Option<String>,
    pub local_start_0based: usize,
    pub local_end_0based_exclusive: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    pub strand: String,
    pub matched_sequence: String,
    pub seed_class: MirnaSeedClass,
    pub region_context_sequence: String,
    pub evidence_tags: Vec<String>,
    pub notes: Vec<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MirnaTargetHitGroup {
    pub region_class: MirnaRegionClass,
    pub transcript_id: Option<String>,
    pub seed_class: MirnaSeedClass,
    pub hits: Vec<MirnaTargetHit>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MirnaTargetScanSummaryCount {
    pub region_class: MirnaRegionClass,
    pub seed_class: MirnaSeedClass,
    pub count: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MirnaTargetScanReport {
    pub schema: String,
    pub query_mirna_id: String,
    pub mature_sequence: String,
    pub accession: Option<String>,
    pub source: String,
    pub target_sequence_id: String,
    pub target_gene_or_sequence_id: String,
    pub seed_motif_table: Vec<MirnaSeedMotif>,
    pub grouped_hits: Vec<MirnaTargetHitGroup>,
    pub summary_counts: Vec<MirnaTargetScanSummaryCount>,
    pub warnings: Vec<String>,
}

#[derive(Clone, Debug)]
struct TranscriptPartition {
    transcript_id: Option<String>,
    feature_id: Option<String>,
    parent_feature: Option<String>,
    gene: Option<String>,
    reverse: bool,
    exons: Vec<(usize, usize)>,
    coding_ranges: Vec<(usize, usize)>,
}

#[derive(Clone, Debug)]
struct ScanRegion {
    class: MirnaRegionClass,
    transcript_id: Option<String>,
    feature_id: Option<String>,
    parent_feature: Option<String>,
    reverse: bool,
    range: (usize, usize),
}

pub fn default_seed_classes() -> Vec<MirnaSeedClass> {
    vec![
        MirnaSeedClass::EightMer,
        MirnaSeedClass::SevenMerM8,
        MirnaSeedClass::SevenMerA1,
        MirnaSeedClass::SixMer,
    ]
}

pub fn default_scan_regions() -> Vec<MirnaRegionClass> {
    vec![
        MirnaRegionClass::ThreePrimeUtr,
        MirnaRegionClass::CodingExon,
        MirnaRegionClass::NoncodingExon,
        MirnaRegionClass::Intron,
        MirnaRegionClass::ExonIntronBoundary,
        MirnaRegionClass::IntronExonBoundary,
    ]
}

pub fn catalog_record(query: &str) -> Option<MirnaCatalogRecord> {
    let normalized = normalize_catalog_key(query);
    let aliases = [
        "hsa-mir-96-5p",
        "hsa-mir-96",
        "mir-96-5p",
        "mimat0000095",
        "407053",
        "ncbi gene 407053",
    ];
    if aliases
        .iter()
        .any(|alias| normalize_catalog_key(alias) == normalized)
    {
        Some(MirnaCatalogRecord {
            schema: MIRNA_CATALOG_RECORD_SCHEMA.to_string(),
            id: "hsa-miR-96-5p".to_string(),
            aliases: vec![
                "MIMAT0000095".to_string(),
                "miR-96-5p".to_string(),
                "NCBI Gene 407053".to_string(),
            ],
            mature_sequence: "UUUGGCACUAGCACAUUUUUGCU".to_string(),
            accession: Some("MIMAT0000095".to_string()),
            source: "GENtle built-in microRNA seed catalog".to_string(),
            ncbi_gene_id: Some("407053".to_string()),
            notes: vec![
                "Rat Tp73 experimental context can be recorded as orthologous evidence for PMID 37099528; human TP73 seed matches remain sequence candidates unless direct human validation is supplied.".to_string(),
            ],
        })
    } else {
        None
    }
}

fn catalog_record_for_mature_sequence(normalized_dna_sequence: &str) -> Option<MirnaCatalogRecord> {
    // The v1 catalog has one built-in entry; iterate catalog rows here when it grows.
    let record = catalog_record("hsa-miR-96-5p")?;
    let record_sequence = normalize_mature_sequence(&record.mature_sequence).ok()?;
    (record_sequence == normalized_dna_sequence).then_some(record)
}

pub fn resolve_mirna_query(
    mirna: &str,
    mature_sequence: Option<&str>,
) -> Result<MirnaQueryRecord, String> {
    if let Some(sequence) = mature_sequence {
        let normalized = normalize_mature_sequence(sequence)?;
        if let Some(record) = catalog_record(mirna)
            .filter(|record| {
                normalize_mature_sequence(&record.mature_sequence)
                    .ok()
                    .as_deref()
                    == Some(normalized.as_str())
            })
            .or_else(|| catalog_record_for_mature_sequence(&normalized))
        {
            return Ok(MirnaQueryRecord {
                id: record.id,
                mature_sequence: normalized,
                accession: record.accession,
                source: format!("direct_sequence_input; matched {}", record.source),
                aliases: record.aliases,
            });
        }
        return Ok(MirnaQueryRecord {
            id: mirna.trim().to_string(),
            mature_sequence: normalized,
            accession: catalog_record(mirna).and_then(|record| record.accession),
            source: "direct_sequence_input".to_string(),
            aliases: catalog_record(mirna)
                .map(|record| record.aliases)
                .unwrap_or_default(),
        });
    }
    if looks_like_sequence(mirna) {
        let normalized = normalize_mature_sequence(mirna)?;
        if let Some(record) = catalog_record_for_mature_sequence(&normalized) {
            return Ok(MirnaQueryRecord {
                id: record.id,
                mature_sequence: normalized,
                accession: record.accession,
                source: format!("direct_sequence_input; matched {}", record.source),
                aliases: record.aliases,
            });
        }
        return Ok(MirnaQueryRecord {
            id: "direct_microRNA_sequence".to_string(),
            mature_sequence: normalized,
            accession: None,
            source: "direct_sequence_input".to_string(),
            aliases: vec![],
        });
    }
    if let Some(record) = catalog_record(mirna) {
        return Ok(MirnaQueryRecord {
            id: record.id,
            mature_sequence: record.mature_sequence,
            accession: record.accession,
            source: record.source,
            aliases: record.aliases,
        });
    }
    Err(format!(
        "Unknown microRNA '{mirna}'. Provide a mature sequence with --mature-sequence."
    ))
}

pub fn explain_seed_motifs(
    mirna: &str,
    mature_sequence: Option<&str>,
) -> Result<serde_json::Value, String> {
    let query = resolve_mirna_query(mirna, mature_sequence)?;
    let motifs = seed_motifs_for_query(&query, &default_seed_classes())?;
    Ok(serde_json::json!({
        "schema": MIRNA_SEED_EXPLANATION_SCHEMA,
        "query_mirna_id": query.id,
        "mature_sequence": query.mature_sequence,
        "accession": query.accession,
        "source": query.source,
        "seed_motif_table": motifs
    }))
}

pub fn scan_mirna_target_sequence(
    sequence: &DNAsequence,
    target_sequence_id: &str,
    target_gene_filter: Option<&str>,
    request: &MirnaTargetScanRequest,
) -> Result<MirnaTargetScanReport, String> {
    let query = resolve_mirna_query(&request.mirna, request.mature_sequence.as_deref())?;
    let seed_motif_table = seed_motifs_for_query(&query, &request.seed_classes)?;
    let requested_regions: BTreeSet<MirnaRegionClass> = request.regions.iter().copied().collect();
    let mut warnings = vec![
        "Seed matches are reported as sequence evidence only; GENtle does not infer functional repression from a seed match alone.".to_string(),
    ];
    let transcripts = transcript_partitions(sequence, target_gene_filter, request, &mut warnings);
    if transcripts.is_empty() {
        warnings.push(format!(
            "No annotated mRNA transcript matched target '{}'.",
            target_gene_filter.unwrap_or(target_sequence_id)
        ));
    }
    let regions = scan_regions_for_transcripts(
        sequence.len(),
        &transcripts,
        &requested_regions,
        request.boundary_flank_bp,
        &mut warnings,
    );
    let evidence_context = evidence_context_notes(&query, target_gene_filter, request);
    let mut flat_hits = vec![];
    for region in regions {
        let Some(region_sequence) =
            oriented_region_sequence(sequence, region.range, region.reverse)
        else {
            continue;
        };
        for motif in &seed_motif_table {
            for offset in find_all(&region_sequence, motif.target_motif.as_bytes()) {
                let (start, end) = map_oriented_hit_to_sequence(
                    region.range,
                    offset,
                    motif.target_motif.len(),
                    region.reverse,
                );
                let mut evidence_tags = vec!["exact_seed_candidate".to_string()];
                let mut notes = request.evidence_notes.clone();
                if let Some(note) = &request.species_note {
                    notes.push(note.clone());
                }
                if let Some(note) = &evidence_context {
                    evidence_tags.push("orthologous_experimental_context".to_string());
                    notes.push(note.clone());
                }
                notes.push(
                    "No direct experimental validation is asserted by this report; results are sequence-evidence candidates only.".to_string(),
                );
                let context = oriented_context_sequence(
                    sequence,
                    region.range,
                    region.reverse,
                    offset,
                    motif.target_motif.len(),
                    DEFAULT_CONTEXT_BP,
                );
                flat_hits.push(MirnaTargetHit {
                    region_class: region.class,
                    transcript_id: region.transcript_id.clone(),
                    feature_id: region.feature_id.clone(),
                    parent_feature: region.parent_feature.clone(),
                    local_start_0based: start,
                    local_end_0based_exclusive: end,
                    genomic_start_1based: start.saturating_add(1),
                    genomic_end_1based: end,
                    strand: if region.reverse {
                        "-".to_string()
                    } else {
                        "+".to_string()
                    },
                    matched_sequence: bytes_to_string(
                        &region_sequence[offset..offset + motif.target_motif.len()],
                    ),
                    seed_class: motif.seed_class,
                    region_context_sequence: context,
                    evidence_tags,
                    notes,
                });
            }
        }
    }
    flat_hits.sort_by(|a, b| {
        a.region_class
            .cmp(&b.region_class)
            .then(a.transcript_id.cmp(&b.transcript_id))
            .then(a.seed_class.cmp(&b.seed_class))
            .then(a.local_start_0based.cmp(&b.local_start_0based))
            .then(
                a.local_end_0based_exclusive
                    .cmp(&b.local_end_0based_exclusive),
            )
    });
    let grouped_hits = grouped_hits(flat_hits);
    let summary_counts = summary_counts(&grouped_hits, &requested_regions, &request.seed_classes);
    Ok(MirnaTargetScanReport {
        schema: MIRNA_TARGET_SCAN_SCHEMA.to_string(),
        query_mirna_id: query.id,
        mature_sequence: query.mature_sequence,
        accession: query.accession,
        source: query.source,
        target_sequence_id: target_sequence_id.to_string(),
        target_gene_or_sequence_id: target_gene_filter.unwrap_or(target_sequence_id).to_string(),
        seed_motif_table,
        grouped_hits,
        summary_counts,
        warnings,
    })
}

fn seed_motifs_for_query(
    query: &MirnaQueryRecord,
    classes: &[MirnaSeedClass],
) -> Result<Vec<MirnaSeedMotif>, String> {
    let mature = normalize_mature_sequence(&query.mature_sequence)?;
    if mature.len() < 8 {
        return Err("microRNA mature sequence must contain at least 8 bases".to_string());
    }
    classes
        .iter()
        .map(|class| {
            let (positions, seed, motif) = match class {
                MirnaSeedClass::EightMer => {
                    let seed = &mature[1..8];
                    ("2..8", seed.to_string(), format!("{}A", revcomp(seed)))
                }
                MirnaSeedClass::SevenMerM8 => {
                    let seed = &mature[1..8];
                    ("2..8", seed.to_string(), revcomp(seed))
                }
                MirnaSeedClass::SevenMerA1 => {
                    let seed = &mature[1..7];
                    ("2..7", seed.to_string(), format!("{}A", revcomp(seed)))
                }
                MirnaSeedClass::SixMer => {
                    let seed = &mature[1..7];
                    ("2..7", seed.to_string(), revcomp(seed))
                }
            };
            Ok(MirnaSeedMotif {
                seed_class: *class,
                mirna_positions_1based: positions.to_string(),
                mirna_seed_5p: seed,
                target_motif: motif,
            })
        })
        .collect()
}

fn transcript_partitions(
    sequence: &DNAsequence,
    target_gene_filter: Option<&str>,
    request: &MirnaTargetScanRequest,
    warnings: &mut Vec<String>,
) -> Vec<TranscriptPartition> {
    let target = target_gene_filter.map(|value| value.to_ascii_lowercase());
    let transcript_filter = request
        .transcript_filter
        .as_deref()
        .map(|value| value.to_ascii_lowercase());
    let features = sequence.features();
    let mut partitions = vec![];
    for (feature_index, feature) in features.iter().enumerate() {
        if !feature.kind.eq_ignore_ascii_case("mRNA") {
            continue;
        }
        let gene = qualifier_first(feature, "gene");
        let transcript_id = qualifier_first(feature, "transcript_id");
        let product = qualifier_first(feature, "product");
        if let Some(target) = &target
            && !qualifier_matches(&gene, target)
            && !qualifier_matches(&transcript_id, target)
            && !qualifier_matches(&product, target)
        {
            continue;
        }
        if let Some(filter) = &transcript_filter {
            let haystacks = [
                transcript_id.as_deref().unwrap_or_default(),
                product.as_deref().unwrap_or_default(),
            ];
            if !haystacks
                .iter()
                .any(|value| value.to_ascii_lowercase().contains(filter))
            {
                continue;
            }
        }
        let mut exons = ranges_for_feature(feature, sequence.len());
        if exons.is_empty() {
            warnings.push(format!(
                "Transcript {} has no parseable exon ranges.",
                transcript_id.as_deref().unwrap_or("<unnamed>")
            ));
            continue;
        }
        exons = merge_ranges(exons);
        let cds_ranges = matching_cds_ranges(features, gene.as_deref(), &exons, sequence.len());
        if cds_ranges.is_empty() {
            warnings.push(format!(
                "Transcript {} has no matching CDS feature; exons are treated as noncoding.",
                transcript_id.as_deref().unwrap_or("<unnamed>")
            ));
        }
        partitions.push(TranscriptPartition {
            transcript_id,
            feature_id: Some(format!("feature:{feature_index}")),
            parent_feature: gene.clone(),
            gene,
            reverse: feature_is_reverse(feature),
            exons,
            coding_ranges: cds_ranges,
        });
    }
    partitions
}

fn scan_regions_for_transcripts(
    seq_len: usize,
    transcripts: &[TranscriptPartition],
    requested: &BTreeSet<MirnaRegionClass>,
    flank: usize,
    warnings: &mut Vec<String>,
) -> Vec<ScanRegion> {
    let mut regions = vec![];
    for transcript in transcripts {
        let coding = merge_ranges(transcript.coding_ranges.clone());
        let three_prime_utr = three_prime_utr_ranges(transcript, &coding);
        if requested.contains(&MirnaRegionClass::ThreePrimeUtr) {
            push_ranges(
                &mut regions,
                transcript,
                MirnaRegionClass::ThreePrimeUtr,
                &three_prime_utr,
            );
        }
        if requested.contains(&MirnaRegionClass::CodingExon) {
            let coding_exons = intersect_ranges(&transcript.exons, &coding);
            push_ranges(
                &mut regions,
                transcript,
                MirnaRegionClass::CodingExon,
                &coding_exons,
            );
        }
        if requested.contains(&MirnaRegionClass::NoncodingExon) {
            let noncoding = subtract_ranges(&transcript.exons, &coding);
            let noncoding_without_3utr = subtract_ranges(&noncoding, &three_prime_utr);
            push_ranges(
                &mut regions,
                transcript,
                MirnaRegionClass::NoncodingExon,
                &noncoding_without_3utr,
            );
        }
        if requested.contains(&MirnaRegionClass::Intron)
            || requested.contains(&MirnaRegionClass::ExonIntronBoundary)
            || requested.contains(&MirnaRegionClass::IntronExonBoundary)
        {
            let introns = intron_ranges(&transcript.exons);
            if introns.is_empty() && transcript.exons.len() > 1 {
                warnings.push(format!(
                    "Transcript {} has adjacent exons with no intronic gap after normalization.",
                    transcript.transcript_id.as_deref().unwrap_or("<unnamed>")
                ));
            }
            if requested.contains(&MirnaRegionClass::Intron) {
                push_ranges(&mut regions, transcript, MirnaRegionClass::Intron, &introns);
            }
            if requested.contains(&MirnaRegionClass::ExonIntronBoundary)
                || requested.contains(&MirnaRegionClass::IntronExonBoundary)
            {
                push_boundary_ranges(&mut regions, transcript, seq_len, flank, requested);
            }
        }
        if requested.contains(&MirnaRegionClass::WholeTranscript) {
            push_ranges(
                &mut regions,
                transcript,
                MirnaRegionClass::WholeTranscript,
                &transcript.exons,
            );
        }
        if requested.contains(&MirnaRegionClass::WholeGene)
            && let (Some(start), Some(end)) = (
                transcript.exons.iter().map(|range| range.0).min(),
                transcript.exons.iter().map(|range| range.1).max(),
            )
        {
            push_ranges(
                &mut regions,
                transcript,
                MirnaRegionClass::WholeGene,
                &[(start, end)],
            );
        }
    }
    regions
}

fn push_ranges(
    regions: &mut Vec<ScanRegion>,
    transcript: &TranscriptPartition,
    class: MirnaRegionClass,
    ranges: &[(usize, usize)],
) {
    for range in ranges.iter().copied().filter(|(start, end)| start < end) {
        regions.push(ScanRegion {
            class,
            transcript_id: transcript.transcript_id.clone(),
            feature_id: transcript.feature_id.clone(),
            parent_feature: transcript
                .parent_feature
                .clone()
                .or_else(|| transcript.gene.clone()),
            reverse: transcript.reverse,
            range,
        });
    }
}

fn push_boundary_ranges(
    regions: &mut Vec<ScanRegion>,
    transcript: &TranscriptPartition,
    seq_len: usize,
    flank: usize,
    requested: &BTreeSet<MirnaRegionClass>,
) {
    let mut exons = transcript.exons.clone();
    exons.sort_unstable();
    for pair in exons.windows(2) {
        let left = pair[0];
        let right = pair[1];
        if left.1 >= right.0 {
            continue;
        }
        if requested.contains(&MirnaRegionClass::ExonIntronBoundary) {
            let range = (
                left.1.saturating_sub(flank),
                left.1.saturating_add(flank).min(seq_len),
            );
            push_ranges(
                regions,
                transcript,
                MirnaRegionClass::ExonIntronBoundary,
                &[range],
            );
        }
        if requested.contains(&MirnaRegionClass::IntronExonBoundary) {
            let range = (
                right.0.saturating_sub(flank),
                right.0.saturating_add(flank).min(seq_len),
            );
            push_ranges(
                regions,
                transcript,
                MirnaRegionClass::IntronExonBoundary,
                &[range],
            );
        }
    }
}

fn grouped_hits(hits: Vec<MirnaTargetHit>) -> Vec<MirnaTargetHitGroup> {
    let mut groups: BTreeMap<
        (MirnaRegionClass, Option<String>, MirnaSeedClass),
        Vec<MirnaTargetHit>,
    > = BTreeMap::new();
    for hit in hits {
        groups
            .entry((hit.region_class, hit.transcript_id.clone(), hit.seed_class))
            .or_default()
            .push(hit);
    }
    groups
        .into_iter()
        .map(
            |((region_class, transcript_id, seed_class), hits)| MirnaTargetHitGroup {
                region_class,
                transcript_id,
                seed_class,
                hits,
            },
        )
        .collect()
}

fn summary_counts(
    groups: &[MirnaTargetHitGroup],
    requested_regions: &BTreeSet<MirnaRegionClass>,
    seed_classes: &[MirnaSeedClass],
) -> Vec<MirnaTargetScanSummaryCount> {
    let mut counts: BTreeMap<(MirnaRegionClass, MirnaSeedClass), usize> = BTreeMap::new();
    for &region_class in requested_regions {
        for &seed_class in seed_classes {
            counts.entry((region_class, seed_class)).or_default();
        }
    }
    for group in groups {
        *counts
            .entry((group.region_class, group.seed_class))
            .or_default() += group.hits.len();
    }
    counts
        .into_iter()
        .map(
            |((region_class, seed_class), count)| MirnaTargetScanSummaryCount {
                region_class,
                seed_class,
                count,
            },
        )
        .collect()
}

fn matching_cds_ranges(
    features: &[Feature],
    gene: Option<&str>,
    transcript_exons: &[(usize, usize)],
    seq_len: usize,
) -> Vec<(usize, usize)> {
    let mut ranges = vec![];
    for feature in features {
        if !feature.kind.eq_ignore_ascii_case("CDS") {
            continue;
        }
        if let Some(gene) = gene
            && !qualifier_values(feature, "gene").any(|value| value.eq_ignore_ascii_case(gene))
        {
            continue;
        }
        let feature_ranges = ranges_for_feature(feature, seq_len);
        if !intersect_ranges(&feature_ranges, transcript_exons).is_empty() {
            ranges.extend(feature_ranges);
        }
    }
    merge_ranges(ranges)
}

fn three_prime_utr_ranges(
    transcript: &TranscriptPartition,
    coding_ranges: &[(usize, usize)],
) -> Vec<(usize, usize)> {
    if coding_ranges.is_empty() {
        return vec![];
    }
    let coding_boundary = if transcript.reverse {
        coding_ranges.iter().map(|range| range.0).min()
    } else {
        coding_ranges.iter().map(|range| range.1).max()
    };
    let Some(boundary) = coding_boundary else {
        return vec![];
    };
    let mut utr = vec![];
    for &(start, end) in &transcript.exons {
        if transcript.reverse {
            if end <= boundary {
                utr.push((start, end));
            } else if start < boundary {
                utr.push((start, boundary));
            }
        } else if start >= boundary {
            utr.push((start, end));
        } else if end > boundary {
            utr.push((boundary, end));
        }
    }
    merge_ranges(utr)
}

fn intron_ranges(exons: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut sorted = exons.to_vec();
    sorted.sort_unstable();
    sorted
        .windows(2)
        .filter_map(|pair| {
            let intron = (pair[0].1, pair[1].0);
            (intron.0 < intron.1).then_some(intron)
        })
        .collect()
}

fn ranges_for_feature(feature: &Feature, seq_len: usize) -> Vec<(usize, usize)> {
    let mut ranges = vec![];
    collect_location_ranges_usize(&feature.location, &mut ranges);
    ranges
        .into_iter()
        .map(|(start, end)| (start.min(seq_len), end.min(seq_len)))
        .filter(|(start, end)| start < end)
        .collect()
}

fn merge_ranges(mut ranges: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    ranges.sort_unstable();
    let mut merged: Vec<(usize, usize)> = vec![];
    for (start, end) in ranges {
        if let Some(last) = merged.last_mut()
            && start <= last.1
        {
            last.1 = last.1.max(end);
            continue;
        }
        if start < end {
            merged.push((start, end));
        }
    }
    merged
}

fn intersect_ranges(left: &[(usize, usize)], right: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut out = vec![];
    for &(a_start, a_end) in left {
        for &(b_start, b_end) in right {
            let start = a_start.max(b_start);
            let end = a_end.min(b_end);
            if start < end {
                out.push((start, end));
            }
        }
    }
    merge_ranges(out)
}

fn subtract_ranges(base: &[(usize, usize)], subtract: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut out = vec![];
    for &(start, end) in base {
        let mut cursor = start;
        for &(sub_start, sub_end) in subtract {
            if sub_end <= cursor || sub_start >= end {
                continue;
            }
            if cursor < sub_start {
                out.push((cursor, sub_start.min(end)));
            }
            cursor = cursor.max(sub_end);
            if cursor >= end {
                break;
            }
        }
        if cursor < end {
            out.push((cursor, end));
        }
    }
    merge_ranges(out)
}

fn oriented_region_sequence(
    sequence: &DNAsequence,
    range: (usize, usize),
    reverse: bool,
) -> Option<Vec<u8>> {
    let bytes = sequence.forward_bytes().get(range.0..range.1)?.to_vec();
    Some(if reverse {
        revcomp_bytes(&bytes)
    } else {
        bytes
    })
}

fn oriented_context_sequence(
    sequence: &DNAsequence,
    range: (usize, usize),
    reverse: bool,
    offset: usize,
    motif_len: usize,
    context_bp: usize,
) -> String {
    let Some(region_sequence) = oriented_region_sequence(sequence, range, reverse) else {
        return String::new();
    };
    let start = offset.saturating_sub(context_bp);
    let end = offset
        .saturating_add(motif_len)
        .saturating_add(context_bp)
        .min(region_sequence.len());
    bytes_to_string(&region_sequence[start..end])
}

fn map_oriented_hit_to_sequence(
    range: (usize, usize),
    offset: usize,
    len: usize,
    reverse: bool,
) -> (usize, usize) {
    if reverse {
        let end = range.1.saturating_sub(offset);
        (end.saturating_sub(len), end)
    } else {
        let start = range.0.saturating_add(offset);
        (start, start.saturating_add(len))
    }
}

fn find_all(haystack: &[u8], needle: &[u8]) -> Vec<usize> {
    if needle.is_empty() || haystack.len() < needle.len() {
        return vec![];
    }
    haystack
        .windows(needle.len())
        .enumerate()
        .filter_map(|(idx, window)| window.eq_ignore_ascii_case(needle).then_some(idx))
        .collect()
}

fn evidence_context_notes(
    query: &MirnaQueryRecord,
    target_gene_filter: Option<&str>,
    request: &MirnaTargetScanRequest,
) -> Option<String> {
    if request
        .evidence_notes
        .iter()
        .any(|note| note.to_ascii_lowercase().contains("pmid 37099528"))
    {
        return Some(
            "Orthologous experimental context supplied by user note: PMID 37099528.".to_string(),
        );
    }
    let target_is_tp73 = target_gene_filter
        .map(|target| target.eq_ignore_ascii_case("TP73") || target.eq_ignore_ascii_case("Tp73"))
        .unwrap_or(false);
    if target_is_tp73 && query.id.eq_ignore_ascii_case("hsa-miR-96-5p") {
        return Some(
            "Orthologous experimental context: rat Tp73 / miR-96 evidence is described in PMID 37099528; this is not direct human TP73 validation."
                .to_string(),
        );
    }
    None
}

fn qualifier_first(feature: &Feature, key: &str) -> Option<String> {
    qualifier_values(feature, key).next().map(ToOwned::to_owned)
}

fn qualifier_values<'a>(feature: &'a Feature, key: &'a str) -> impl Iterator<Item = &'a str> {
    feature.qualifier_values(key).map(str::trim)
}

fn qualifier_matches(value: &Option<String>, target: &str) -> bool {
    value
        .as_deref()
        .map(|value| value.to_ascii_lowercase().contains(target))
        .unwrap_or(false)
}

fn looks_like_sequence(raw: &str) -> bool {
    let trimmed = raw.trim();
    !trimmed.is_empty()
        && trimmed
            .chars()
            .all(|c| matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T' | 'U'))
}

fn normalize_mature_sequence(raw: &str) -> Result<String, String> {
    let normalized: String = raw
        .chars()
        .filter(|c| !c.is_ascii_whitespace())
        .map(|c| match c.to_ascii_uppercase() {
            'U' => 'T',
            other => other,
        })
        .collect();
    if normalized.is_empty() {
        return Err("microRNA mature sequence is empty".to_string());
    }
    if normalized
        .chars()
        .any(|c| !matches!(c, 'A' | 'C' | 'G' | 'T'))
    {
        return Err(format!(
            "microRNA mature sequence contains unsupported characters: '{raw}'"
        ));
    }
    Ok(normalized)
}

fn normalize_catalog_key(raw: &str) -> String {
    raw.trim()
        .to_ascii_lowercase()
        .replace(['_', '-'], "")
        .replace(' ', "")
}

fn revcomp(sequence: &str) -> String {
    bytes_to_string(&revcomp_bytes(sequence.as_bytes()))
}

fn revcomp_bytes(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            _ => b'N',
        })
        .collect()
}

fn bytes_to_string(bytes: &[u8]) -> String {
    bytes
        .iter()
        .map(|base| base.to_ascii_uppercase() as char)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna_sequence;

    #[test]
    fn hsa_mir_96_seed_motifs_match_expected_target_sites() {
        let query = resolve_mirna_query("hsa-miR-96-5p", None).expect("query");
        let motifs = seed_motifs_for_query(&query, &default_seed_classes()).expect("motifs");
        let by_class: BTreeMap<MirnaSeedClass, String> = motifs
            .into_iter()
            .map(|motif| (motif.seed_class, motif.target_motif))
            .collect();
        assert_eq!(by_class[&MirnaSeedClass::EightMer], "GTGCCAAA");
        assert_eq!(by_class[&MirnaSeedClass::SevenMerM8], "GTGCCAA");
        assert_eq!(by_class[&MirnaSeedClass::SevenMerA1], "TGCCAAA");
        assert_eq!(by_class[&MirnaSeedClass::SixMer], "TGCCAA");
    }

    #[test]
    fn direct_mature_sequence_resolves_matching_catalog_record() {
        let query = resolve_mirna_query("custom-label", Some("UUUGGCACUAGCACAUUUUUGCU"))
            .expect("catalog match");
        assert_eq!(query.id, "hsa-miR-96-5p");
        assert_eq!(query.accession.as_deref(), Some("MIMAT0000095"));
        let raw_query = resolve_mirna_query("UUUGGCACUAGCACAUUUUUGCU", None).expect("raw match");
        assert_eq!(raw_query.id, "hsa-miR-96-5p");
    }

    #[test]
    fn helper_outputs_emit_declared_schemas() {
        let explanation = explain_seed_motifs("hsa-miR-96-5p", None).expect("explain seed");
        assert_eq!(
            explanation.get("schema").and_then(|value| value.as_str()),
            Some(MIRNA_SEED_EXPLANATION_SCHEMA)
        );
        let record = catalog_record("hsa-miR-96-5p").expect("catalog record");
        assert_eq!(record.schema, MIRNA_CATALOG_RECORD_SCHEMA);
    }

    #[test]
    fn generic_target_scan_notes_do_not_leak_tp73_wording() {
        let mut sequence =
            DNAsequence::from_sequence("ATGAAACCCCTGCCAAATTT").expect("synthetic sequence");
        sequence.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: gb_io::seq::Location::simple_range(0, 20),
            qualifiers: vec![
                ("gene".into(), Some("EGFR".to_string())),
                ("transcript_id".into(), Some("EGFR_SYNTHETIC.1".to_string())),
            ],
        });
        sequence.features_mut().push(gb_io::seq::Feature {
            kind: "CDS".into(),
            location: gb_io::seq::Location::simple_range(0, 9),
            qualifiers: vec![("gene".into(), Some("EGFR".to_string()))],
        });
        let request = MirnaTargetScanRequest::with_defaults("hsa-miR-96-5p");
        let report =
            scan_mirna_target_sequence(&sequence, "egfr_synthetic", Some("EGFR"), &request)
                .expect("scan");
        let notes = report
            .grouped_hits
            .iter()
            .flat_map(|group| group.hits.iter())
            .flat_map(|hit| hit.notes.iter())
            .cloned()
            .collect::<Vec<_>>();
        assert!(
            notes
                .iter()
                .any(|note| note.contains("No direct experimental validation"))
        );
        assert!(
            notes.iter().all(|note| !note.contains("TP73")),
            "generic target notes must not mention TP73: {notes:?}"
        );
    }

    #[test]
    fn tp73_refseq_scan_reports_3utr_mir_96_candidate_without_validation_claim() {
        let sequence = dna_sequence::load_from_file("test_files/tp73.ncbi.gb").expect("load TP73");
        let request = MirnaTargetScanRequest::with_defaults("hsa-miR-96-5p");
        let report = scan_mirna_target_sequence(&sequence, "tp73.ncbi.gb", Some("TP73"), &request)
            .expect("scan");
        assert_eq!(report.schema, MIRNA_TARGET_SCAN_SCHEMA);
        assert!(
            report
                .warnings
                .iter()
                .any(|warning| warning.contains("does not infer functional repression"))
        );
        assert!(
            report
                .grouped_hits
                .iter()
                .flat_map(|group| group.hits.iter())
                .any(|hit| hit.region_class == MirnaRegionClass::ThreePrimeUtr
                    && hit.seed_class == MirnaSeedClass::SevenMerA1
                    && hit.matched_sequence == "TGCCAAA"
                    && hit.local_start_0based == 80976
                    && hit.local_end_0based_exclusive == 80983)
        );
        assert!(
            report
                .grouped_hits
                .iter()
                .any(|group| group.region_class == MirnaRegionClass::Intron)
        );
        assert!(
            report
                .summary_counts
                .iter()
                .any(|count| count.region_class == MirnaRegionClass::CodingExon)
        );
        assert!(report.summary_counts.iter().any(|count| count.region_class
            == MirnaRegionClass::ExonIntronBoundary
            || count.region_class == MirnaRegionClass::IntronExonBoundary));
        assert!(report.summary_counts.iter().any(|count| count.region_class
            == MirnaRegionClass::ThreePrimeUtr
            && count.seed_class == MirnaSeedClass::EightMer
            && count.count == 0));
        assert!(
            report
                .grouped_hits
                .iter()
                .flat_map(|group| group.hits.iter())
                .flat_map(|hit| hit.evidence_tags.iter())
                .all(|tag| tag != "directly_validated_human_site")
        );
        assert!(
            report
                .grouped_hits
                .iter()
                .flat_map(|group| group.hits.iter())
                .any(|hit| hit
                    .evidence_tags
                    .iter()
                    .any(|tag| tag == "orthologous_experimental_context")
                    && hit.notes.iter().any(|note| note.contains("PMID 37099528")))
        );
    }
}
