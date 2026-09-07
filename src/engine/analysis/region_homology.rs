//! Read-only local-BLAST homology screens for portable genomic regions.

use super::*;
use gentle_protocol as gp;
use std::sync::{LazyLock, Mutex};

static REGION_HOMOLOGY_CACHE: LazyLock<
    Mutex<BTreeMap<String, gp::GenomicRegionHomologyScreenReport>>,
> = LazyLock::new(|| Mutex::new(BTreeMap::new()));

fn homology_error(code: ErrorCode, message: impl Into<String>) -> EngineError {
    EngineError::new(code, message)
}

fn canonical_digest<T: Serialize>(value: &T, label: &str) -> Result<String, EngineError> {
    let json = serde_json::to_string(value).map_err(|error| {
        homology_error(
            ErrorCode::Internal,
            format!("could not serialize {label} for deterministic digest: {error}"),
        )
    })?;
    Ok(sha256_prefixed_str(&json))
}

fn region_homology_cache_key(
    query_sequence_sha256: &str,
    effective: &gp::GenomicRegionHomologyEffectiveRequest,
    database_cache_identity: &str,
) -> Result<String, EngineError> {
    Ok(sha256_prefixed_str(&format!(
        "{}\0{}\0{}\0{}",
        query_sequence_sha256,
        canonical_digest(effective, "effective request")?,
        gp::GENOMIC_REGION_HOMOLOGY_PROJECTION_VERSION,
        database_cache_identity
    )))
}

pub(crate) fn validate_genomic_region_homology_report(
    report: &gp::GenomicRegionHomologyScreenReport,
) -> Result<(), EngineError> {
    if report.schema != gp::GENOMIC_REGION_HOMOLOGY_SCREEN_SCHEMA {
        return Err(homology_error(
            ErrorCode::InvalidInput,
            format!("unsupported homology report schema '{}'", report.schema),
        ));
    }
    if report.projection_version != gp::GENOMIC_REGION_HOMOLOGY_PROJECTION_VERSION {
        return Err(homology_error(
            ErrorCode::InvalidInput,
            format!(
                "unsupported homology projection version '{}'",
                report.projection_version
            ),
        ));
    }
    let mut content = report.clone();
    let expected = content.content_sha256.clone();
    content.content_sha256.clear();
    content.op_id = None;
    content.run_id = None;
    if expected.is_empty() || canonical_digest(&content, "homology report")? != expected {
        return Err(homology_error(
            ErrorCode::InvalidInput,
            "homology report content digest does not match its payload",
        ));
    }
    Ok(())
}

fn reference_contig_aliases(reference: &gp::GenomicRegionReference) -> BTreeSet<String> {
    std::iter::once(reference.contig_name.as_str())
        .chain(reference.contig_accession.as_deref())
        .chain(reference.contig_aliases.iter().map(String::as_str))
        .flat_map(|value| {
            let mut aliases = vec![value.trim().to_string()];
            let trimmed = value.trim();
            if let Some(rest) = trimmed.strip_prefix("chr") {
                aliases.push(rest.to_string());
            } else if !trimmed.is_empty() {
                aliases.push(format!("chr{trimmed}"));
            }
            aliases
        })
        .filter(|value| !value.is_empty())
        .collect()
}

fn subject_alias_map(
    catalog: &GenomeCatalog,
    genome_id: &str,
    cache_dir: Option<&str>,
) -> BTreeMap<String, String> {
    let mut aliases = BTreeMap::new();
    let records = catalog
        .list_chromosome_lengths(genome_id, cache_dir)
        .unwrap_or_default();
    let mut versionless = BTreeMap::<String, Vec<String>>::new();
    for record in records {
        let canonical = record.chromosome;
        aliases.insert(canonical.clone(), canonical.clone());
        aliases.insert(canonical.to_ascii_lowercase(), canonical.clone());
        if let Some((base, version)) = canonical.rsplit_once('.')
            && version.chars().all(|ch| ch.is_ascii_digit())
        {
            versionless
                .entry(base.to_string())
                .or_default()
                .push(canonical.clone());
        }
    }
    for (base, values) in versionless {
        if values.len() == 1 {
            aliases.insert(base.clone(), values[0].clone());
            aliases.insert(base.to_ascii_lowercase(), values[0].clone());
        }
    }
    aliases
}

fn normalize_subject_id(raw: &str, aliases: &BTreeMap<String, String>) -> String {
    let raw = raw.trim();
    let mut candidates = vec![raw];
    if let Some(first) = raw.split_whitespace().next() {
        candidates.push(first);
    }
    candidates.extend(
        raw.split('|')
            .map(str::trim)
            .filter(|value| !value.is_empty()),
    );
    for candidate in candidates {
        if let Some(canonical) = aliases
            .get(candidate)
            .or_else(|| aliases.get(&candidate.to_ascii_lowercase()))
        {
            return canonical.clone();
        }
    }
    raw.to_string()
}

fn references_same_assembly(
    reference: &gp::GenomicRegionReference,
    inspection: &BlastDatabaseInspectionReport,
) -> bool {
    inspection.source_genome_id == reference.assembly_name
        || inspection
            .source_assembly
            .as_deref()
            .is_some_and(|assembly| {
                assembly.eq_ignore_ascii_case(&reference.assembly_name)
                    || reference
                        .assembly_accession
                        .as_deref()
                        .is_some_and(|accession| assembly.eq_ignore_ascii_case(accession))
            })
}

fn database_binding(
    inspection: &BlastDatabaseInspectionReport,
) -> gp::GenomicRegionHomologyDatabaseBinding {
    gp::GenomicRegionHomologyDatabaseBinding {
        genome_id: inspection.source_genome_id.clone(),
        index_kind: inspection.index_kind.as_str().to_string(),
        source_assembly: inspection.source_assembly.clone(),
        source_release: inspection.source_release.clone(),
        masking: inspection.masking.clone(),
        prefix: inspection.prefix.clone(),
        blast_database_version: inspection.blast_database_version.clone(),
        sequence_count: inspection.sequence_count,
        total_bases: inspection.total_bases,
        tool_executable: inspection.tool_executable.clone(),
        tool_version: inspection.tool_version.clone(),
        content_fingerprint: inspection.content_fingerprint.clone(),
        fingerprint_algorithm: inspection.fingerprint_algorithm.clone(),
        validation_status: inspection.validation_status.clone(),
    }
}

fn references_same_contig(reference: &gp::GenomicRegionReference, subject_id: &str) -> bool {
    let subject_tokens = std::iter::once(subject_id)
        .chain(subject_id.split('|'))
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .flat_map(|value| {
            let mut values = vec![value.to_ascii_lowercase()];
            if let Some(rest) = value.strip_prefix("chr") {
                values.push(rest.to_ascii_lowercase());
            } else {
                values.push(format!("chr{}", value.to_ascii_lowercase()));
            }
            values
        })
        .collect::<BTreeSet<_>>();
    reference_contig_aliases(reference)
        .into_iter()
        .map(|value| value.to_ascii_lowercase())
        .any(|value| subject_tokens.contains(&value))
}

fn intervals_overlap(left_start: u64, left_end: u64, right_start: u64, right_end: u64) -> bool {
    left_start < right_end && right_start < left_end
}

fn classify_locus(
    query: &gp::GenomicRegionHomologyQueryBinding,
    target: &gp::GenomicRegionHomologyTargetRequest,
    subject_id: &str,
    subject_start: u64,
    subject_end: u64,
) -> (gp::GenomicRegionHomologyLocusClass, Option<String>) {
    if target.role == gp::GenomicRegionHomologyTargetRole::SameGenome {
        let interval = &query.region.interval;
        if references_same_contig(&interval.reference, subject_id)
            && intervals_overlap(
                interval.start_0based,
                interval.end_0based_exclusive,
                subject_start,
                subject_end,
            )
        {
            return (gp::GenomicRegionHomologyLocusClass::SameGenomeSelf, None);
        }
        return (gp::GenomicRegionHomologyLocusClass::SameGenomeNonself, None);
    }
    for expected in &target.expected_loci {
        if references_same_contig(&expected.reference, subject_id)
            && intervals_overlap(
                expected.start_0based,
                expected.end_0based_exclusive,
                subject_start,
                subject_end,
            )
        {
            return (
                gp::GenomicRegionHomologyLocusClass::ExpectedOrtholog,
                Some(expected.evidence_id.clone()),
            );
        }
    }
    (
        gp::GenomicRegionHomologyLocusClass::CrossSpeciesUnassigned,
        None,
    )
}

fn convert_blast_hit(
    target_genome_id: &str,
    hit: &crate::genomes::BlastHit,
    aliases: &BTreeMap<String, String>,
    ordinal: usize,
) -> Option<gp::GenomicRegionHomologyHsp> {
    let aligned_query = hit.aligned_query.as_ref()?.to_ascii_uppercase();
    let aligned_subject = hit.aligned_subject.as_ref()?.to_ascii_uppercase();
    if aligned_query.len() != aligned_subject.len() {
        return None;
    }
    let subject_id = normalize_subject_id(&hit.subject_id, aliases);
    let strand = if hit.subject_start <= hit.subject_end {
        gp::GenomicRegionStrand::Plus
    } else {
        gp::GenomicRegionStrand::Minus
    };
    let query_start_0based = hit.query_start.min(hit.query_end).saturating_sub(1);
    let query_end_0based_exclusive = hit.query_start.max(hit.query_end);
    let subject_start_0based = hit.subject_start.min(hit.subject_end).saturating_sub(1) as u64;
    let subject_end_0based_exclusive = hit.subject_start.max(hit.subject_end) as u64;
    let identity = format!(
        "{}\0{}\0{}\0{}\0{}\0{}\0{}\0{}",
        target_genome_id,
        subject_id,
        query_start_0based,
        query_end_0based_exclusive,
        subject_start_0based,
        subject_end_0based_exclusive,
        hit.bit_score.to_bits(),
        ordinal
    );
    Some(gp::GenomicRegionHomologyHsp {
        hsp_id: short_sha256_id("region_hsp", &identity),
        target_genome_id: target_genome_id.to_string(),
        subject_id_raw: hit.subject_id.clone(),
        subject_id,
        strand,
        query_start_0based,
        query_end_0based_exclusive,
        subject_start_0based,
        subject_end_0based_exclusive,
        identity_percent: hit.identity_percent,
        alignment_length_bp: hit.alignment_length,
        mismatches: hit.mismatches,
        gap_opens: hit.gap_opens,
        evalue: hit.evalue,
        bit_score: hit.bit_score,
        aligned_query,
        aligned_subject,
    })
}

fn oriented_subject_start(hsp: &gp::GenomicRegionHomologyHsp) -> i128 {
    match hsp.strand {
        gp::GenomicRegionStrand::Minus => -(hsp.subject_end_0based_exclusive as i128),
        _ => hsp.subject_start_0based as i128,
    }
}

fn oriented_subject_interval(hsp: &gp::GenomicRegionHomologyHsp) -> (i128, i128) {
    match hsp.strand {
        gp::GenomicRegionStrand::Minus => (
            -(hsp.subject_end_0based_exclusive as i128),
            -(hsp.subject_start_0based as i128),
        ),
        _ => (
            hsp.subject_start_0based as i128,
            hsp.subject_end_0based_exclusive as i128,
        ),
    }
}

fn hsp_chain_gap(
    left: &gp::GenomicRegionHomologyHsp,
    right: &gp::GenomicRegionHomologyHsp,
) -> Option<usize> {
    if left.subject_id != right.subject_id || left.strand != right.strand {
        return None;
    }
    let (left_subject_start, left_subject_end) = oriented_subject_interval(left);
    let (right_subject_start, _) = oriented_subject_interval(right);
    if right.query_start_0based < left.query_start_0based
        || right_subject_start < left_subject_start
    {
        return None;
    }
    let query_overlaps = right.query_start_0based < left.query_end_0based_exclusive;
    let subject_overlaps = right_subject_start < left_subject_end;
    if query_overlaps != subject_overlaps {
        return None;
    }
    if query_overlaps {
        return Some(0);
    }
    let query_gap = right.query_start_0based - left.query_end_0based_exclusive;
    let subject_gap = right_subject_start.saturating_sub(left_subject_end) as u128;
    usize::try_from(subject_gap.max(query_gap as u128)).ok()
}

fn chain_hsps(
    hsps: &[gp::GenomicRegionHomologyHsp],
    max_gap_bp: usize,
) -> Vec<Vec<gp::GenomicRegionHomologyHsp>> {
    let mut sorted = hsps.to_vec();
    sorted.sort_by(|left, right| {
        (
            left.subject_id.as_str(),
            left.strand.bed_value(),
            left.query_start_0based,
            oriented_subject_start(left),
            left.hsp_id.as_str(),
        )
            .cmp(&(
                right.subject_id.as_str(),
                right.strand.bed_value(),
                right.query_start_0based,
                oriented_subject_start(right),
                right.hsp_id.as_str(),
            ))
    });
    let mut chains: Vec<Vec<gp::GenomicRegionHomologyHsp>> = vec![];
    for hsp in sorted {
        let append_to = chains
            .iter()
            .enumerate()
            .filter_map(|(index, chain)| {
                hsp_chain_gap(chain.last()?, &hsp)
                    .filter(|gap| *gap <= max_gap_bp)
                    .map(|gap| (gap, index))
            })
            .min();
        if let Some((_, index)) = append_to {
            chains[index].push(hsp);
        } else {
            chains.push(vec![hsp]);
        }
    }
    chains
}

fn locus_from_chain(
    query: &gp::GenomicRegionHomologyQueryBinding,
    target: &gp::GenomicRegionHomologyTargetRequest,
    chain: &[gp::GenomicRegionHomologyHsp],
) -> gp::GenomicRegionHomologyLocus {
    let subject_id = chain[0].subject_id.clone();
    let strand = chain[0].strand;
    let query_start = chain
        .iter()
        .map(|hsp| hsp.query_start_0based)
        .min()
        .unwrap_or(0);
    let query_end = chain
        .iter()
        .map(|hsp| hsp.query_end_0based_exclusive)
        .max()
        .unwrap_or(query_start);
    let subject_start = chain
        .iter()
        .map(|hsp| hsp.subject_start_0based)
        .min()
        .unwrap_or(0);
    let subject_end = chain
        .iter()
        .map(|hsp| hsp.subject_end_0based_exclusive)
        .max()
        .unwrap_or(subject_start);
    let (locus_class, orthology_evidence_id) =
        classify_locus(query, target, &subject_id, subject_start, subject_end);
    let aligned_total = chain
        .iter()
        .map(|hsp| hsp.alignment_length_bp)
        .sum::<usize>()
        .max(1);
    let identity_percent = chain
        .iter()
        .map(|hsp| hsp.identity_percent * hsp.alignment_length_bp as f64)
        .sum::<f64>()
        / aligned_total as f64;
    let bit_score = chain.iter().map(|hsp| hsp.bit_score).sum::<f64>();
    let identity = format!(
        "{}\0{}\0{}\0{}\0{}\0{}\0{}",
        target.genome_id,
        subject_id,
        strand.bed_value(),
        query_start,
        query_end,
        subject_start,
        subject_end
    );
    gp::GenomicRegionHomologyLocus {
        locus_id: short_sha256_id("region_locus", &identity),
        target_genome_id: target.genome_id.clone(),
        target_role: target.role,
        locus_class,
        subject_id,
        strand,
        query_start_0based: query_start,
        query_end_0based_exclusive: query_end,
        subject_start_0based: subject_start,
        subject_end_0based_exclusive: subject_end,
        identity_percent,
        query_coverage_percent: (query_end.saturating_sub(query_start) as f64
            / query.sequence.len().max(1) as f64)
            * 100.0,
        bit_score,
        source_hsp_ids: chain.iter().map(|hsp| hsp.hsp_id.clone()).collect(),
        orthology_evidence_id,
    }
}

fn locus_class_rank(class: gp::GenomicRegionHomologyLocusClass) -> u8 {
    match class {
        gp::GenomicRegionHomologyLocusClass::Query => 0,
        gp::GenomicRegionHomologyLocusClass::ExpectedOrtholog => 1,
        gp::GenomicRegionHomologyLocusClass::CrossSpeciesUnassigned => 2,
        gp::GenomicRegionHomologyLocusClass::SameGenomeSelf => 3,
        gp::GenomicRegionHomologyLocusClass::SameGenomeNonself => 4,
    }
}

fn project_locus_row(
    query_sequence: &str,
    locus: &gp::GenomicRegionHomologyLocus,
    hsps: &BTreeMap<String, gp::GenomicRegionHomologyHsp>,
) -> (
    gp::GenomicRegionHomologyAlignmentRow,
    Vec<gp::GenomicRegionHomologyOmittedInsertion>,
) {
    #[derive(Clone)]
    struct ProjectedBase {
        symbol: u8,
        hsp_id: String,
    }
    let mut projection: Vec<Option<ProjectedBase>> = vec![None; query_sequence.len()];
    let mut conflicts = vec![];
    let mut insertions = vec![];
    let mut ordered_hsps = locus
        .source_hsp_ids
        .iter()
        .filter_map(|id| hsps.get(id))
        .collect::<Vec<_>>();
    ordered_hsps.sort_by(|left, right| {
        right
            .bit_score
            .total_cmp(&left.bit_score)
            .then(right.alignment_length_bp.cmp(&left.alignment_length_bp))
            .then(left.evalue.total_cmp(&right.evalue))
            .then(left.hsp_id.cmp(&right.hsp_id))
    });
    for hsp in ordered_hsps {
        let mut query_position = hsp.query_start_0based;
        let mut subject_position = match hsp.strand {
            gp::GenomicRegionStrand::Minus => hsp.subject_end_0based_exclusive,
            _ => hsp.subject_start_0based,
        };
        let query_bytes = hsp.aligned_query.as_bytes();
        let subject_bytes = hsp.aligned_subject.as_bytes();
        let mut column = 0usize;
        while column < query_bytes.len() {
            if query_bytes[column] == b'-' {
                let insertion_start_column = column;
                let target_start = subject_position;
                let mut sequence = String::new();
                while column < query_bytes.len() && query_bytes[column] == b'-' {
                    if subject_bytes[column] != b'-' {
                        sequence.push(subject_bytes[column] as char);
                        match hsp.strand {
                            gp::GenomicRegionStrand::Minus => {
                                subject_position = subject_position.saturating_sub(1)
                            }
                            _ => subject_position = subject_position.saturating_add(1),
                        }
                    }
                    column += 1;
                }
                let target_end = subject_position;
                let (target_start_0based, target_end_0based_exclusive) =
                    if hsp.strand == gp::GenomicRegionStrand::Minus {
                        (target_end, target_start)
                    } else {
                        (target_start, target_end)
                    };
                let identity = format!(
                    "{}\0{}\0{}\0{}\0{}",
                    locus.locus_id, hsp.hsp_id, query_position, insertion_start_column, sequence
                );
                insertions.push(gp::GenomicRegionHomologyOmittedInsertion {
                    insertion_id: short_sha256_id("region_insertion", &identity),
                    alignment_row_id: locus.locus_id.clone(),
                    hsp_id: hsp.hsp_id.clone(),
                    query_anchor_0based: query_position.min(query_sequence.len()),
                    target_start_0based,
                    target_end_0based_exclusive,
                    length_bp: sequence.len(),
                    target_sequence: sequence,
                    strand: hsp.strand,
                });
                continue;
            }
            if query_position >= projection.len() {
                break;
            }
            let query_base = query_bytes[column].to_ascii_uppercase();
            let subject_base = subject_bytes[column].to_ascii_uppercase();
            let symbol = if subject_base == b'-' {
                b'-'
            } else if query_base == subject_base {
                b'.'
            } else {
                subject_base
            };
            let projected = ProjectedBase {
                symbol,
                hsp_id: hsp.hsp_id.clone(),
            };
            if let Some(retained) = &projection[query_position] {
                if retained.symbol != symbol {
                    conflicts.push(gp::GenomicRegionHomologyAlignmentConflict {
                        query_position_0based: query_position,
                        retained_hsp_id: retained.hsp_id.clone(),
                        rejected_hsp_id: hsp.hsp_id.clone(),
                        retained_symbol: (retained.symbol as char).to_string(),
                        rejected_symbol: (symbol as char).to_string(),
                    });
                }
            } else {
                projection[query_position] = Some(projected);
            }
            query_position += 1;
            if subject_base != b'-' {
                match hsp.strand {
                    gp::GenomicRegionStrand::Minus => {
                        subject_position = subject_position.saturating_sub(1)
                    }
                    _ => subject_position = subject_position.saturating_add(1),
                }
            }
            column += 1;
        }
    }
    let query_projection = projection
        .iter()
        .map(|base| base.as_ref().map_or(' ', |base| base.symbol as char))
        .collect::<String>();
    let exact_match_count = projection
        .iter()
        .filter(|base| base.as_ref().is_some_and(|base| base.symbol == b'.'))
        .count();
    let covered_query_base_count = projection.iter().filter(|base| base.is_some()).count();
    let row_id = short_sha256_id(
        "region_alignment",
        &format!("{}\0{}", locus.locus_id, query_projection),
    );
    for insertion in &mut insertions {
        insertion.alignment_row_id = row_id.clone();
    }
    (
        gp::GenomicRegionHomologyAlignmentRow {
            row_id,
            locus_id: locus.locus_id.clone(),
            target_genome_id: locus.target_genome_id.clone(),
            locus_class: locus.locus_class,
            subject_id: locus.subject_id.clone(),
            strand: locus.strand,
            subject_start_0based: locus.subject_start_0based,
            subject_end_0based_exclusive: locus.subject_end_0based_exclusive,
            query_projection,
            exact_match_count,
            covered_query_base_count,
            omitted_insertion_ids: insertions
                .iter()
                .map(|insertion| insertion.insertion_id.clone())
                .collect(),
            source_hsp_ids: locus.source_hsp_ids.clone(),
            conflicts,
        },
        insertions,
    )
}

fn block_genomic_interval(
    query: &gp::GenomicRegionHomologyQueryBinding,
    start: usize,
    end: usize,
) -> (u64, u64) {
    let interval = &query.region.interval;
    if interval.strand == gp::GenomicRegionStrand::Minus {
        (
            interval.end_0based_exclusive.saturating_sub(end as u64),
            interval.end_0based_exclusive.saturating_sub(start as u64),
        )
    } else {
        (
            interval.start_0based.saturating_add(start as u64),
            interval.start_0based.saturating_add(end as u64),
        )
    }
}

fn support_class_for_row(
    row: &gp::GenomicRegionHomologyAlignmentRow,
) -> Option<gp::GenomicRegionHomologySupportClass> {
    match row.locus_class {
        gp::GenomicRegionHomologyLocusClass::ExpectedOrtholog => {
            Some(gp::GenomicRegionHomologySupportClass::ExpectedOrtholog)
        }
        gp::GenomicRegionHomologyLocusClass::CrossSpeciesUnassigned => {
            Some(gp::GenomicRegionHomologySupportClass::CrossSpeciesUnassigned)
        }
        gp::GenomicRegionHomologyLocusClass::SameGenomeNonself => {
            Some(gp::GenomicRegionHomologySupportClass::SameGenomeNonself)
        }
        _ => None,
    }
}

fn target_support_class(
    target: &gp::GenomicRegionHomologyTargetRequest,
) -> gp::GenomicRegionHomologySupportClass {
    match target.role {
        gp::GenomicRegionHomologyTargetRole::ExpectedOrtholog => {
            gp::GenomicRegionHomologySupportClass::ExpectedOrtholog
        }
        gp::GenomicRegionHomologyTargetRole::SameGenome => {
            gp::GenomicRegionHomologySupportClass::SameGenomeNonself
        }
        gp::GenomicRegionHomologyTargetRole::CrossSpeciesUnassigned => {
            gp::GenomicRegionHomologySupportClass::CrossSpeciesUnassigned
        }
    }
}

fn call_conserved_blocks(
    query: &gp::GenomicRegionHomologyQueryBinding,
    effective: &gp::GenomicRegionHomologyEffectiveRequest,
    targets: &[gp::GenomicRegionHomologyTargetResult],
    rows: &[gp::GenomicRegionHomologyAlignmentRow],
) -> Vec<gp::GenomicRegionHomologyConservedBlock> {
    let hsp_ids_by_row = rows
        .iter()
        .map(|row| (row.row_id.as_str(), row.source_hsp_ids.as_slice()))
        .collect::<BTreeMap<_, _>>();
    let mut blocks = vec![];
    for class in [
        gp::GenomicRegionHomologySupportClass::ExpectedOrtholog,
        gp::GenomicRegionHomologySupportClass::CrossSpeciesUnassigned,
        gp::GenomicRegionHomologySupportClass::SameGenomeNonself,
    ] {
        let available_genome_ids = targets
            .iter()
            .filter(|target| target_support_class(&target.target) == class)
            .filter(|target| {
                matches!(
                    target.status,
                    gp::GenomicRegionHomologyTargetStatus::Available
                        | gp::GenomicRegionHomologyTargetStatus::NoAcceptedSimilarity
                )
            })
            .map(|target| target.target.genome_id.clone())
            .collect::<BTreeSet<_>>();
        let unavailable_genome_ids = targets
            .iter()
            .filter(|target| target_support_class(&target.target) == class)
            .filter(|target| {
                matches!(
                    target.status,
                    gp::GenomicRegionHomologyTargetStatus::Unavailable
                        | gp::GenomicRegionHomologyTargetStatus::SearchOutputTooBroad
                )
            })
            .map(|target| target.target.genome_id.clone())
            .collect::<Vec<_>>();
        let class_rows = rows
            .iter()
            .filter(|row| support_class_for_row(row) == Some(class))
            .collect::<Vec<_>>();
        let mut start = 0usize;
        while start < query.sequence.len() {
            let supporters = class_rows
                .iter()
                .filter(|row| row.query_projection.as_bytes().get(start) == Some(&b'.'))
                .map(|row| row.row_id.clone())
                .collect::<BTreeSet<_>>();
            if supporters.is_empty() {
                start += 1;
                continue;
            }
            let mut end = start + 1;
            while end < query.sequence.len() {
                let next = class_rows
                    .iter()
                    .filter(|row| row.query_projection.as_bytes().get(end) == Some(&b'.'))
                    .map(|row| row.row_id.clone())
                    .collect::<BTreeSet<_>>();
                if next != supporters {
                    break;
                }
                end += 1;
            }
            if end - start >= effective.policy.min_conserved_block_bp {
                let supporting_genome_ids = class_rows
                    .iter()
                    .filter(|row| supporters.contains(&row.row_id))
                    .map(|row| row.target_genome_id.clone())
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect::<Vec<_>>();
                let source_hsp_ids = supporters
                    .iter()
                    .flat_map(|row_id| {
                        hsp_ids_by_row
                            .get(row_id.as_str())
                            .into_iter()
                            .flat_map(|ids| ids.iter().cloned())
                    })
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect::<Vec<_>>();
                let (genomic_start, genomic_end) = block_genomic_interval(query, start, end);
                let identity = format!(
                    "{}\0{}\0{}\0{}\0{}",
                    query.region.content_sha256,
                    class.as_str(),
                    start,
                    end,
                    supporters.iter().cloned().collect::<Vec<_>>().join(",")
                );
                blocks.push(gp::GenomicRegionHomologyConservedBlock {
                    block_id: short_sha256_id("conserved_block", &identity),
                    support_class: class,
                    query_start_0based: start,
                    query_end_0based_exclusive: end,
                    query_sequence: query.sequence[start..end].to_string(),
                    genomic_start_0based: genomic_start,
                    genomic_end_0based_exclusive: genomic_end,
                    genomic_strand: query.region.interval.strand,
                    supporting_row_ids: supporters.into_iter().collect(),
                    supporting_genome_ids: supporting_genome_ids.clone(),
                    available_genome_ids: available_genome_ids.iter().cloned().collect(),
                    unavailable_genome_ids: unavailable_genome_ids.clone(),
                    support_fraction: if available_genome_ids.is_empty() {
                        0.0
                    } else {
                        supporting_genome_ids.len() as f64 / available_genome_ids.len() as f64
                    },
                    mean_identity_percent: 100.0,
                    source_hsp_ids,
                });
            }
            start = end;
        }
    }
    blocks.sort_by(|left, right| {
        (
            left.support_class,
            left.query_start_0based,
            left.query_end_0based_exclusive,
            left.block_id.as_str(),
        )
            .cmp(&(
                right.support_class,
                right.query_start_0based,
                right.query_end_0based_exclusive,
                right.block_id.as_str(),
            ))
    });
    blocks
}

fn finalize_projection(
    query: gp::GenomicRegionHomologyQueryBinding,
    effective_request: gp::GenomicRegionHomologyEffectiveRequest,
    mut targets: Vec<gp::GenomicRegionHomologyTargetResult>,
    mut hsps: Vec<gp::GenomicRegionHomologyHsp>,
    request_sha256: String,
) -> Result<gp::GenomicRegionHomologyScreenReport, EngineError> {
    targets.sort_by(|left, right| left.target.genome_id.cmp(&right.target.genome_id));
    hsps.sort_by(|left, right| {
        (
            left.target_genome_id.as_str(),
            left.subject_id.as_str(),
            left.query_start_0based,
            left.subject_start_0based,
            left.hsp_id.as_str(),
        )
            .cmp(&(
                right.target_genome_id.as_str(),
                right.subject_id.as_str(),
                right.query_start_0based,
                right.subject_start_0based,
                right.hsp_id.as_str(),
            ))
    });
    let hsp_lookup = hsps
        .iter()
        .map(|hsp| (hsp.hsp_id.clone(), hsp.clone()))
        .collect::<BTreeMap<_, _>>();
    let target_lookup = effective_request
        .targets
        .iter()
        .map(|target| (target.genome_id.as_str(), target))
        .collect::<BTreeMap<_, _>>();
    let mut loci = vec![];
    for (target_id, target) in &target_lookup {
        let target_hsps = hsps
            .iter()
            .filter(|hsp| hsp.target_genome_id == **target_id)
            .cloned()
            .collect::<Vec<_>>();
        let mut target_loci = chain_hsps(&target_hsps, effective_request.policy.max_chain_gap_bp)
            .into_iter()
            .filter(|chain| !chain.is_empty())
            .map(|chain| locus_from_chain(&query, target, &chain))
            .collect::<Vec<_>>();
        target_loci.sort_by(|left, right| {
            locus_class_rank(left.locus_class)
                .cmp(&locus_class_rank(right.locus_class))
                .then(
                    right
                        .query_coverage_percent
                        .total_cmp(&left.query_coverage_percent),
                )
                .then(right.bit_score.total_cmp(&left.bit_score))
                .then(left.subject_id.cmp(&right.subject_id))
                .then(left.subject_start_0based.cmp(&right.subject_start_0based))
                .then(left.locus_id.cmp(&right.locus_id))
        });
        target_loci.truncate(effective_request.policy.max_loci_per_target);
        if let Some(result) = targets
            .iter_mut()
            .find(|result| result.target.genome_id == **target_id)
        {
            result.retained_locus_count = target_loci.len();
            if result.status == gp::GenomicRegionHomologyTargetStatus::Available
                && target_loci.is_empty()
            {
                result.status = gp::GenomicRegionHomologyTargetStatus::NoAcceptedSimilarity;
            }
        }
        loci.extend(target_loci);
    }
    loci.sort_by(|left, right| {
        locus_class_rank(left.locus_class)
            .cmp(&locus_class_rank(right.locus_class))
            .then(left.target_genome_id.cmp(&right.target_genome_id))
            .then(left.subject_id.cmp(&right.subject_id))
            .then(left.subject_start_0based.cmp(&right.subject_start_0based))
            .then(left.locus_id.cmp(&right.locus_id))
    });
    let query_row = gp::GenomicRegionHomologyAlignmentRow {
        row_id: short_sha256_id("region_alignment_query", &query.sequence_sha256),
        locus_id: "query".to_string(),
        target_genome_id: effective_request.query_genome_id.clone(),
        locus_class: gp::GenomicRegionHomologyLocusClass::Query,
        subject_id: query.region.interval.reference.contig_name.clone(),
        strand: query.region.interval.strand,
        subject_start_0based: query.region.interval.start_0based,
        subject_end_0based_exclusive: query.region.interval.end_0based_exclusive,
        query_projection: query.sequence.clone(),
        exact_match_count: query.sequence.len(),
        covered_query_base_count: query.sequence.len(),
        omitted_insertion_ids: vec![],
        source_hsp_ids: vec![],
        conflicts: vec![],
    };
    let mut rows = vec![query_row];
    let mut insertions = vec![];
    for locus in &loci {
        let (row, mut row_insertions) = project_locus_row(&query.sequence, locus, &hsp_lookup);
        rows.push(row);
        insertions.append(&mut row_insertions);
    }
    let conserved_blocks = call_conserved_blocks(&query, &effective_request, &targets, &rows);
    let mut same_genome_coverage = vec![false; query.sequence.len()];
    for row in rows
        .iter()
        .filter(|row| row.locus_class == gp::GenomicRegionHomologyLocusClass::SameGenomeNonself)
    {
        for (position, symbol) in row.query_projection.bytes().enumerate() {
            if symbol != b' ' {
                same_genome_coverage[position] = true;
            }
        }
    }
    let same_genome_nonself_locus_count = loci
        .iter()
        .filter(|locus| locus.locus_class == gp::GenomicRegionHomologyLocusClass::SameGenomeNonself)
        .count();
    let same_genome_nonself_query_coverage_percent = if query.sequence.is_empty() {
        0.0
    } else {
        same_genome_coverage.iter().filter(|value| **value).count() as f64
            / query.sequence.len() as f64
            * 100.0
    };
    let effective_request_sha256 = canonical_digest(&effective_request, "effective request")?;
    let no_targets = effective_request.targets.is_empty();
    let mut report = gp::GenomicRegionHomologyScreenReport {
        schema: gp::GENOMIC_REGION_HOMOLOGY_SCREEN_SCHEMA.to_string(),
        projection_version: gp::GENOMIC_REGION_HOMOLOGY_PROJECTION_VERSION.to_string(),
        request_sha256,
        effective_request_sha256,
        content_sha256: String::new(),
        query,
        effective_request,
        targets,
        hsps,
        loci,
        alignment_rows: rows,
        omitted_insertions: insertions,
        conserved_blocks,
        same_genome_nonself_locus_count,
        same_genome_nonself_query_coverage_percent,
        warnings: if no_targets {
            vec![
                "No validated local genomic-DNA BLAST indexes were available; the report contains the bound query only."
                    .to_string(),
            ]
        } else {
            vec![]
        },
        non_claims: vec![
            "BLAST similarity alone is not an orthology assertion.".to_string(),
            "Conservation does not prove that a sequence fragment is independently functional."
                .to_string(),
            "Same-genome similarity is reported as repetition or paralog-like similarity, not as a functional off-target claim."
                .to_string(),
            "An unavailable genome is not counted as biological absence of conservation."
                .to_string(),
        ],
        op_id: None,
        run_id: None,
    };
    let mut content = report.clone();
    content.content_sha256.clear();
    content.op_id = None;
    content.run_id = None;
    report.content_sha256 = canonical_digest(&content, "homology report")?;
    Ok(report)
}

fn resolve_query_binding(
    engine: &GentleEngine,
    request: &gp::GenomicRegionHomologyScreenRequest,
    catalog: &GenomeCatalog,
) -> Result<(gp::GenomicRegionHomologyQueryBinding, String), EngineError> {
    if request.set_id.trim().is_empty() || request.region_id.trim().is_empty() {
        return Err(homology_error(
            ErrorCode::InvalidInput,
            "set_id and region_id are required",
        ));
    }
    let store = engine.genomic_region_store_snapshot()?;
    let set = store
        .sets
        .iter()
        .find(|set| set.set_id == request.set_id)
        .ok_or_else(|| {
            homology_error(
                ErrorCode::NotFound,
                format!("genomic region set '{}' not found", request.set_id),
            )
        })?;
    let region = set
        .regions
        .iter()
        .find(|region| region.region_id == request.region_id)
        .cloned()
        .ok_or_else(|| {
            homology_error(
                ErrorCode::NotFound,
                format!(
                    "genomic region '{}' not found in set '{}'",
                    request.region_id, request.set_id
                ),
            )
        })?;
    if request
        .expected_region_content_sha256
        .as_deref()
        .is_some_and(|expected| expected != region.content_sha256)
    {
        return Err(homology_error(
            ErrorCode::InvalidInput,
            "saved genomic region changed after the homology request was prepared",
        ));
    }
    let cache_dir = request
        .cache_dir
        .as_deref()
        .map(str::trim)
        .filter(|v| !v.is_empty());
    if let Some(projection) = region.local_projection.as_ref()
        && projection.status == gp::GenomicRegionLocalProjectionStatus::Current
        && let Some(dna) = engine.state.sequences.get(&projection.seq_id)
        && sha256_prefixed_bytes(dna.forward_bytes()) == projection.sequence_sha256
    {
        let start = usize::try_from(projection.local_start_0based).map_err(|_| {
            homology_error(
                ErrorCode::InvalidInput,
                "local ROI start exceeds platform range",
            )
        })?;
        let end = usize::try_from(projection.local_end_0based_exclusive).map_err(|_| {
            homology_error(
                ErrorCode::InvalidInput,
                "local ROI end exceeds platform range",
            )
        })?;
        if start < end && end <= dna.forward_bytes().len() {
            let mut sequence =
                String::from_utf8_lossy(&dna.forward_bytes()[start..end]).to_ascii_uppercase();
            if projection.local_strand == gp::GenomicRegionStrand::Minus {
                sequence = GentleEngine::reverse_complement(&sequence);
            }
            let query_genome_id = request
                .query_genome_id
                .clone()
                .unwrap_or_else(|| projection.source_genome_id.clone());
            let projection_seq_id = projection.seq_id.clone();
            let projection_sequence_sha256 = projection.sequence_sha256.clone();
            return Ok((
                gp::GenomicRegionHomologyQueryBinding {
                    set_id: set.set_id.clone(),
                    region,
                    sequence_sha256: sha256_prefixed_str(&sequence),
                    sequence,
                    source_kind: "current_local_projection".to_string(),
                    source_resource_id: projection_seq_id,
                    source_fingerprint: Some(projection_sequence_sha256),
                },
                query_genome_id,
            ));
        }
    }
    let query_genome_id = request.query_genome_id.as_deref().ok_or_else(|| {
        homology_error(
            ErrorCode::InvalidInput,
            "query_genome_id is required when the saved region has no current local projection",
        )
    })?;
    let start_1based = usize::try_from(region.interval.start_0based.saturating_add(1))
        .map_err(|_| homology_error(ErrorCode::InvalidInput, "ROI start exceeds platform range"))?;
    let end_1based = usize::try_from(region.interval.end_0based_exclusive)
        .map_err(|_| homology_error(ErrorCode::InvalidInput, "ROI end exceeds platform range"))?;
    let mut sequence = catalog
        .get_sequence_region_with_cache(
            query_genome_id,
            &region.interval.reference.contig_name,
            start_1based,
            end_1based,
            cache_dir,
        )
        .map_err(|error| {
            homology_error(
                ErrorCode::InvalidInput,
                format!("could not resolve saved region sequence: {error}"),
            )
        })?
        .to_ascii_uppercase();
    if region.interval.strand == gp::GenomicRegionStrand::Minus {
        sequence = GentleEngine::reverse_complement(&sequence);
    }
    let inspection = catalog
        .inspect_blast_database(query_genome_id, cache_dir)
        .map_err(|error| homology_error(ErrorCode::InvalidInput, error))?;
    Ok((
        gp::GenomicRegionHomologyQueryBinding {
            set_id: set.set_id.clone(),
            region,
            sequence_sha256: sha256_prefixed_str(&sequence),
            sequence,
            source_kind: "prepared_genome".to_string(),
            source_resource_id: query_genome_id.to_string(),
            source_fingerprint: inspection.and_then(|item| item.content_fingerprint),
        },
        query_genome_id.to_string(),
    ))
}

fn validate_policy(policy: &gp::GenomicRegionHomologySearchPolicy) -> Result<(), EngineError> {
    if !(0.0..=100.0).contains(&policy.min_identity_percent)
        || !policy.min_identity_percent.is_finite()
        || policy.min_alignment_length_bp == 0
        || !policy.max_evalue.is_finite()
        || policy.max_evalue < 0.0
        || policy.max_loci_per_target == 0
        || policy.max_hsps_per_target == 0
        || policy.min_conserved_block_bp == 0
    {
        return Err(homology_error(
            ErrorCode::InvalidInput,
            "homology policy requires finite identity/e-value thresholds and positive alignment, locus, HSP, and conserved-block limits",
        ));
    }
    Ok(())
}

fn emit_region_homology_progress(
    on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    request: &gp::GenomicRegionHomologyScreenRequest,
    phase: &str,
    target_genome_id: Option<&str>,
    target_ordinal: usize,
    target_count: usize,
    detail: impl Into<String>,
    done: bool,
) -> Result<(), EngineError> {
    if on_progress(OperationProgress::GenomicRegionHomology(
        GenomicRegionHomologyProgress {
            set_id: request.set_id.clone(),
            region_id: request.region_id.clone(),
            phase: phase.to_string(),
            target_genome_id: target_genome_id.map(str::to_string),
            target_ordinal,
            target_count,
            detail: detail.into(),
            done,
        },
    )) {
        Ok(())
    } else {
        Err(homology_error(
            ErrorCode::Internal,
            "genomic-region homology screen cancelled during progress reporting",
        ))
    }
}

struct HomologyTargetPreflight {
    target: gp::GenomicRegionHomologyTargetRequest,
    inspection: Option<BlastDatabaseInspectionReport>,
    unavailable_reason: Option<String>,
}

impl GentleEngine {
    pub(crate) fn screen_genomic_region_homology(
        &self,
        request: gp::GenomicRegionHomologyScreenRequest,
        op_id: &str,
        run_id: &str,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<gp::GenomicRegionHomologyScreenReport, EngineError> {
        emit_region_homology_progress(
            on_progress,
            &request,
            "preflight",
            None,
            0,
            request.targets.len(),
            "validating the saved query region and local genomic indexes",
            false,
        )?;
        validate_policy(&request.policy)?;
        let request_sha256 = canonical_digest(&request, "homology request")?;
        let (catalog, catalog_origin) =
            Self::open_reference_genome_catalog(request.catalog_path.as_deref())?;
        let (query, query_genome_id) = resolve_query_binding(self, &request, &catalog)?;
        if query.sequence.is_empty() {
            return Err(homology_error(
                ErrorCode::InvalidInput,
                "saved genomic region resolved to an empty query sequence",
            ));
        }
        let cache_dir = request
            .cache_dir
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty());
        let mut inspections = BTreeMap::new();
        let mut effective_targets = if request.targets.is_empty() {
            let mut resolved = vec![];
            for genome_id in catalog.list_genomes() {
                let Ok(Some(inspection)) = catalog.inspect_blast_database(&genome_id, cache_dir)
                else {
                    continue;
                };
                if inspection.validation_status != "valid"
                    || inspection.index_kind != BlastDatabaseIndexKind::GenomicDna
                {
                    continue;
                }
                let role = if genome_id == query_genome_id
                    || references_same_assembly(&query.region.interval.reference, &inspection)
                {
                    gp::GenomicRegionHomologyTargetRole::SameGenome
                } else {
                    gp::GenomicRegionHomologyTargetRole::CrossSpeciesUnassigned
                };
                inspections.insert(genome_id.clone(), inspection);
                resolved.push(gp::GenomicRegionHomologyTargetRequest {
                    genome_id,
                    required: false,
                    role,
                    expected_loci: vec![],
                });
            }
            resolved
        } else {
            request.targets.clone()
        };
        effective_targets.sort_by(|left, right| left.genome_id.cmp(&right.genome_id));
        let mut dedup = BTreeSet::new();
        for target in &effective_targets {
            if target.genome_id.trim().is_empty() || !dedup.insert(target.genome_id.as_str()) {
                return Err(homology_error(
                    ErrorCode::InvalidInput,
                    "homology targets require distinct, non-empty genome_id values",
                ));
            }
        }
        let effective = gp::GenomicRegionHomologyEffectiveRequest {
            set_id: request.set_id.clone(),
            region_id: request.region_id.clone(),
            region_content_sha256: query.region.content_sha256.clone(),
            query_genome_id,
            targets: effective_targets.clone(),
            catalog_origin,
            cache_dir: request.cache_dir.clone(),
            policy: request.policy.clone(),
        };
        let target_count = effective_targets.len();
        let mut preflight_targets = vec![];
        for (target_index, target) in effective_targets.iter().enumerate() {
            let (inspection, inspection_error) =
                if let Some(inspection) = inspections.remove(&target.genome_id) {
                    (Some(inspection), None)
                } else {
                    match catalog.inspect_blast_database(&target.genome_id, cache_dir) {
                        Ok(inspection) => (inspection, None),
                        Err(error) => (None, Some(error)),
                    }
                };
            let usable = inspection.as_ref().is_some_and(|inspection| {
                inspection.validation_status == "valid"
                    && inspection.index_kind == BlastDatabaseIndexKind::GenomicDna
                    && inspection.content_fingerprint.is_some()
            });
            let unavailable_reason = if usable {
                None
            } else {
                Some(match (inspection_error, inspection.as_ref()) {
                    (Some(error), _) => error,
                    (None, None) => "no prepared BLAST index was found".to_string(),
                    (_, Some(item)) if item.index_kind != BlastDatabaseIndexKind::GenomicDna => {
                        format!(
                            "index kind '{}' is not genomic DNA",
                            item.index_kind.as_str()
                        )
                    }
                    (_, Some(item)) => format!(
                        "BLAST index validation status is '{}' or its content fingerprint is unavailable",
                        item.validation_status
                    ),
                })
            };
            if let Some(detail) = unavailable_reason.as_deref() {
                if target.required {
                    return Err(homology_error(
                        ErrorCode::InvalidInput,
                        format!(
                            "required homology target '{}' is unavailable: {detail}",
                            target.genome_id
                        ),
                    ));
                }
            }
            emit_region_homology_progress(
                on_progress,
                &request,
                "preflight",
                Some(&target.genome_id),
                target_index + 1,
                target_count,
                unavailable_reason
                    .as_deref()
                    .map_or("validated genomic-DNA index", |reason| reason),
                false,
            )?;
            preflight_targets.push(HomologyTargetPreflight {
                target: target.clone(),
                inspection,
                unavailable_reason,
            });
        }
        let database_cache_identity = preflight_targets
            .iter()
            .map(|item| match item.inspection.as_ref() {
                Some(database) => format!(
                    "{}:{}:{}:{}:{}",
                    item.target.genome_id,
                    database.content_fingerprint.as_deref().unwrap_or("missing"),
                    database.tool_version.as_deref().unwrap_or("unknown"),
                    database.validation_status,
                    database.index_kind.as_str(),
                ),
                None => format!(
                    "{}:unavailable:{}",
                    item.target.genome_id,
                    item.unavailable_reason.as_deref().unwrap_or("not_found")
                ),
            })
            .collect::<Vec<_>>()
            .join("\0");
        let cache_key = region_homology_cache_key(
            &query.sequence_sha256,
            &effective,
            &database_cache_identity,
        )?;
        if let Ok(cache) = REGION_HOMOLOGY_CACHE.lock()
            && let Some(cached) = cache.get(&cache_key)
        {
            let mut report = cached.clone();
            report.op_id = Some(op_id.to_string());
            report.run_id = Some(run_id.to_string());
            emit_region_homology_progress(
                on_progress,
                &request,
                "complete",
                None,
                target_count,
                target_count,
                "reused the content-bound local homology result cache",
                true,
            )?;
            return Ok(report);
        }

        let mut target_results = vec![];
        let mut all_hsps = vec![];
        for (target_index, preflight) in preflight_targets.into_iter().enumerate() {
            let target = preflight.target;
            if let Some(detail) = preflight.unavailable_reason {
                target_results.push(gp::GenomicRegionHomologyTargetResult {
                    target,
                    status: gp::GenomicRegionHomologyTargetStatus::Unavailable,
                    database: preflight.inspection.as_ref().map(database_binding),
                    warnings: vec![detail],
                    ..Default::default()
                });
                continue;
            }
            let inspection = preflight.inspection.ok_or_else(|| {
                homology_error(
                    ErrorCode::Internal,
                    format!(
                        "homology target '{}' passed availability checks without an inspection report",
                        target.genome_id
                    ),
                )
            })?;
            emit_region_homology_progress(
                on_progress,
                &request,
                "blast",
                Some(&target.genome_id),
                target_index + 1,
                target_count,
                "searching the validated local genomic database",
                false,
            )?;
            let aliases = subject_alias_map(&catalog, &target.genome_id, cache_dir);
            let blast = catalog
                .blast_sequence_complete_aligned_with_cache(
                    &target.genome_id,
                    &query.sequence,
                    effective.policy.max_hsps_per_target.saturating_add(1),
                    Some("blastn"),
                    cache_dir,
                )
                .map_err(|error| {
                    homology_error(
                        ErrorCode::InvalidInput,
                        format!("BLAST failed for target '{}': {error}", target.genome_id),
                    )
                })?;
            if blast.hits.len() > effective.policy.max_hsps_per_target {
                let warning = format!(
                    "target '{}' produced {} HSPs, exceeding the processing budget of {}",
                    target.genome_id,
                    blast.hits.len(),
                    effective.policy.max_hsps_per_target
                );
                if target.required {
                    return Err(homology_error(ErrorCode::InvalidInput, warning));
                }
                target_results.push(gp::GenomicRegionHomologyTargetResult {
                    target: target.clone(),
                    status: gp::GenomicRegionHomologyTargetStatus::SearchOutputTooBroad,
                    database: Some(database_binding(&inspection)),
                    raw_hsp_count: blast.hits.len(),
                    warnings: vec![warning],
                    ..Default::default()
                });
                continue;
            }
            let mut missing_alignment_strings = 0usize;
            let accepted = blast
                .hits
                .iter()
                .enumerate()
                .filter(|(_, hit)| {
                    hit.identity_percent >= effective.policy.min_identity_percent
                        && hit.alignment_length >= effective.policy.min_alignment_length_bp
                        && hit.evalue <= effective.policy.max_evalue
                })
                .filter_map(|(ordinal, hit)| {
                    let converted = convert_blast_hit(&target.genome_id, hit, &aliases, ordinal);
                    if converted.is_none() {
                        missing_alignment_strings += 1;
                    }
                    converted
                })
                .collect::<Vec<_>>();
            let mut warnings = blast.warnings;
            if missing_alignment_strings > 0 {
                warnings.push(format!(
                    "{missing_alignment_strings} accepted legacy HSP(s) lacked aligned strings and could not be projected"
                ));
            }
            target_results.push(gp::GenomicRegionHomologyTargetResult {
                target: target.clone(),
                status: gp::GenomicRegionHomologyTargetStatus::Available,
                database: Some(database_binding(&inspection)),
                raw_hsp_count: blast.hits.len(),
                accepted_hsp_count: accepted.len(),
                retained_locus_count: 0,
                warnings,
            });
            all_hsps.extend(accepted);
            emit_region_homology_progress(
                on_progress,
                &request,
                "blast",
                Some(&target.genome_id),
                target_index + 1,
                target_count,
                "completed local BLAST search and retained accepted HSPs",
                false,
            )?;
        }
        emit_region_homology_progress(
            on_progress,
            &request,
            "projection",
            None,
            target_count,
            target_count,
            "projecting accepted alignments onto the query coordinate system",
            false,
        )?;
        let mut report =
            finalize_projection(query, effective, target_results, all_hsps, request_sha256)?;
        report.op_id = Some(op_id.to_string());
        report.run_id = Some(run_id.to_string());
        if let Ok(mut cache) = REGION_HOMOLOGY_CACHE.lock() {
            let mut cached = report.clone();
            cached.op_id = None;
            cached.run_id = None;
            cache.insert(cache_key, cached);
        }
        emit_region_homology_progress(
            on_progress,
            &request,
            "block_calling",
            None,
            target_count,
            target_count,
            "completed exact-support block calling and report composition",
            true,
        )?;
        Ok(report)
    }

    pub(crate) fn assess_promoter_conserved_modules(
        &self,
        request: gp::PromoterModuleAssessmentRequest,
        op_id: &str,
        run_id: &str,
    ) -> Result<gp::PromoterModuleAssessmentReport, EngineError> {
        let homology = request.homology_report.as_ref();
        validate_genomic_region_homology_report(homology)?;
        if !request.max_same_genome_query_coverage_percent.is_finite()
            || !(0.0..=100.0).contains(&request.max_same_genome_query_coverage_percent)
        {
            return Err(homology_error(
                ErrorCode::InvalidInput,
                "max_same_genome_query_coverage_percent must be between 0 and 100",
            ));
        }
        let mut evidence = request.selected_evidence_spans.clone();
        evidence.sort_by(|left, right| left.evidence_id.cmp(&right.evidence_id));
        let query_len = homology.query.sequence.len();
        let mut ids = BTreeSet::new();
        for span in &evidence {
            if span.evidence_id.trim().is_empty()
                || !ids.insert(span.evidence_id.as_str())
                || span.query_start_0based >= span.query_end_0based_exclusive
                || span.query_end_0based_exclusive > query_len
            {
                return Err(homology_error(
                    ErrorCode::InvalidInput,
                    "selected evidence requires unique ids and non-empty query-bounded intervals",
                ));
            }
        }
        let expected_blocks = homology
            .conserved_blocks
            .iter()
            .filter(|block| {
                block.support_class == gp::GenomicRegionHomologySupportClass::ExpectedOrtholog
            })
            .collect::<Vec<_>>();
        let required = evidence
            .iter()
            .filter(|span| span.required)
            .collect::<Vec<_>>();
        let fully_containing = expected_blocks.iter().find(|block| {
            required.iter().all(|span| {
                block.query_start_0based <= span.query_start_0based
                    && block.query_end_0based_exclusive >= span.query_end_0based_exclusive
            })
        });
        let mut covering_blocks = BTreeSet::new();
        let all_required_covered = required.iter().all(|span| {
            let covering = expected_blocks.iter().find(|block| {
                block.query_start_0based <= span.query_start_0based
                    && block.query_end_0based_exclusive >= span.query_end_0based_exclusive
            });
            if let Some(block) = covering {
                covering_blocks.insert(block.block_id.clone());
                true
            } else {
                false
            }
        });
        let mut selected_blocks = covering_blocks.into_iter().collect::<Vec<_>>();
        selected_blocks.sort_by_key(|id| {
            homology
                .conserved_blocks
                .iter()
                .find(|block| &block.block_id == id)
                .map_or(usize::MAX, |block| block.query_start_0based)
        });
        let spacing_retained = selected_blocks.windows(2).all(|pair| {
            let left = homology
                .conserved_blocks
                .iter()
                .find(|block| block.block_id == pair[0]);
            let right = homology
                .conserved_blocks
                .iter()
                .find(|block| block.block_id == pair[1]);
            matches!((left, right), (Some(left), Some(right)) if right.query_start_0based.saturating_sub(left.query_end_0based_exclusive) <= request.max_partner_gap_bp)
        });
        let ortholog_targets_available = homology.targets.iter().any(|target| {
            target.target.role == gp::GenomicRegionHomologyTargetRole::ExpectedOrtholog
                && matches!(
                    target.status,
                    gp::GenomicRegionHomologyTargetStatus::Available
                        | gp::GenomicRegionHomologyTargetStatus::NoAcceptedSimilarity
                )
        });
        let repetitive = homology.same_genome_nonself_query_coverage_percent
            >= request.max_same_genome_query_coverage_percent;
        let same_genome_targets = homology
            .effective_request
            .targets
            .iter()
            .filter(|target| target.role == gp::GenomicRegionHomologyTargetRole::SameGenome)
            .collect::<Vec<_>>();
        let same_genome_assessed = !same_genome_targets.is_empty()
            && same_genome_targets.iter().all(|requested| {
                homology.targets.iter().any(|result| {
                    result.target == **requested
                        && matches!(
                            result.status,
                            gp::GenomicRegionHomologyTargetStatus::Available
                                | gp::GenomicRegionHomologyTargetStatus::NoAcceptedSimilarity
                        )
                })
            });
        let hypothesis = if !ortholog_targets_available || required.is_empty() {
            gp::PromoterModuleHypothesisKind::InsufficientEvidence
        } else if repetitive {
            gp::PromoterModuleHypothesisKind::RepetitiveOrAmbiguous
        } else if !same_genome_assessed {
            gp::PromoterModuleHypothesisKind::InsufficientEvidence
        } else if fully_containing.is_some() {
            gp::PromoterModuleHypothesisKind::StandaloneReporterCandidate
        } else if all_required_covered && selected_blocks.len() > 1 && spacing_retained {
            gp::PromoterModuleHypothesisKind::PairedContextCandidate
        } else {
            gp::PromoterModuleHypothesisKind::InsufficientEvidence
        };
        if let Some(block) = fully_containing {
            selected_blocks = vec![block.block_id.clone()];
        }
        let decision_trace = vec![
            gp::PromoterModuleDecisionRule {
                rule_id: "ortholog_evidence_available".to_string(),
                description: "At least one explicitly designated ortholog target was evaluable."
                    .to_string(),
                satisfied: ortholog_targets_available,
                detail: format!("{} expected-ortholog block(s)", expected_blocks.len()),
                ..Default::default()
            },
            gp::PromoterModuleDecisionRule {
                rule_id: "selected_evidence_complete".to_string(),
                description: "Every required selected evidence span is represented.".to_string(),
                satisfied: !required.is_empty(),
                evidence_ids: required.iter().map(|span| span.evidence_id.clone()).collect(),
                detail: format!("{} required evidence span(s)", required.len()),
                ..Default::default()
            },
            gp::PromoterModuleDecisionRule {
                rule_id: "single_conserved_block_contains_selected_evidence".to_string(),
                description: "One expected-ortholog exact-support block contains every required evidence span."
                    .to_string(),
                satisfied: fully_containing.is_some(),
                block_ids: fully_containing
                    .iter()
                    .map(|block| block.block_id.clone())
                    .collect(),
                detail: "Containment is geometric and does not prove autonomous function."
                    .to_string(),
                ..Default::default()
            },
            gp::PromoterModuleDecisionRule {
                rule_id: "paired_blocks_preserve_context".to_string(),
                description: "Multiple expected-ortholog blocks jointly cover the evidence while retaining order and allowed spacing."
                    .to_string(),
                satisfied: all_required_covered && selected_blocks.len() > 1 && spacing_retained,
                block_ids: selected_blocks.clone(),
                detail: format!("maximum allowed partner gap {} bp", request.max_partner_gap_bp),
                ..Default::default()
            },
            gp::PromoterModuleDecisionRule {
                rule_id: "same_genome_evidence_available".to_string(),
                description: "Every requested same-genome search completed within its HSP budget."
                    .to_string(),
                satisfied: same_genome_assessed,
                detail: if same_genome_assessed {
                    "assessed under the declared similarity filters".to_string()
                } else {
                    "unassessed: request and complete a same-genome search before interpreting absence of repetition".to_string()
                },
                ..Default::default()
            },
            gp::PromoterModuleDecisionRule {
                rule_id: "same_genome_interpretation_unique".to_string(),
                description: "Same-genome non-self similarity stays below the ambiguity threshold."
                    .to_string(),
                satisfied: same_genome_assessed && !repetitive,
                detail: if same_genome_assessed { format!(
                    "observed {:.3}% versus threshold {:.3}%",
                    homology.same_genome_nonself_query_coverage_percent,
                    request.max_same_genome_query_coverage_percent
                ) } else { "unassessed; a numeric zero is not evidence of uniqueness".to_string() },
                ..Default::default()
            },
        ];
        let alternative_fragments = homology
            .conserved_blocks
            .iter()
            .filter(|block| {
                block.support_class == gp::GenomicRegionHomologySupportClass::ExpectedOrtholog
            })
            .map(|block| gp::PromoterModuleAlternativeFragment {
                fragment_id: format!("fragment_{}", block.block_id),
                query_start_0based: block.query_start_0based,
                query_end_0based_exclusive: block.query_end_0based_exclusive,
                block_ids: vec![block.block_id.clone()],
                covered_evidence_ids: required
                    .iter()
                    .filter(|span| {
                        block.query_start_0based <= span.query_start_0based
                            && block.query_end_0based_exclusive >= span.query_end_0based_exclusive
                    })
                    .map(|span| span.evidence_id.clone())
                    .collect(),
                rationale:
                    "Exact-support block retained as a separately testable fragment hypothesis."
                        .to_string(),
            })
            .collect::<Vec<_>>();
        let request_sha256 = canonical_digest(&request, "promoter-module request")?;
        let mut report = gp::PromoterModuleAssessmentReport {
            schema: gp::PROMOTER_MODULE_ASSESSMENT_SCHEMA.to_string(),
            request_sha256,
            content_sha256: String::new(),
            homology_report_sha256: homology.content_sha256.clone(),
            hypothesis,
            selected_evidence_spans: evidence,
            selected_block_ids: selected_blocks,
            decision_trace,
            alternative_fragments,
            suggested_validation: vec![
                "Compare each candidate block alone with its partner block alone and the combined ordered fragment."
                    .to_string(),
                "Include motif-disrupted controls while preserving fragment length and recording the exact edit."
                    .to_string(),
            ],
            warnings: vec![],
            non_claims: vec![
                "Conservation does not prove autonomous promoter or enhancer activity.".to_string(),
                "A module hypothesis is a traceable reporter-fragment design aid, not a causal regulatory conclusion."
                    .to_string(),
            ],
            op_id: Some(op_id.to_string()),
            run_id: Some(run_id.to_string()),
        };
        let mut content = report.clone();
        content.content_sha256.clear();
        content.op_id = None;
        content.run_id = None;
        report.content_sha256 = canonical_digest(&content, "promoter-module report")?;
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn query() -> gp::GenomicRegionHomologyQueryBinding {
        let sequence = "AACCGGTTAACCGGTT".to_string();
        gp::GenomicRegionHomologyQueryBinding {
            set_id: "synthetic_set".to_string(),
            region: gp::GenomicRegionOfInterest {
                schema: gp::GENOMIC_REGION_OF_INTEREST_SCHEMA.to_string(),
                region_id: "synthetic_region".to_string(),
                interval: gp::GenomicRegionInterval {
                    reference: gp::GenomicRegionReference {
                        species_scientific_name: Some("Synthetic organism".to_string()),
                        taxon_id: Some(1),
                        assembly_name: "synthetic_v1".to_string(),
                        contig_name: "chrQ".to_string(),
                        ..Default::default()
                    },
                    start_0based: 100,
                    end_0based_exclusive: 116,
                    strand: gp::GenomicRegionStrand::Plus,
                    ..Default::default()
                },
                content_sha256: "sha256:query".to_string(),
                identity_sha256: "sha256:identity".to_string(),
                ..Default::default()
            },
            sequence_sha256: sha256_prefixed_str(&sequence),
            sequence,
            source_kind: "synthetic".to_string(),
            source_resource_id: "fixture".to_string(),
            source_fingerprint: Some("sha256:fixture".to_string()),
        }
    }

    fn target(role: gp::GenomicRegionHomologyTargetRole) -> gp::GenomicRegionHomologyTargetRequest {
        gp::GenomicRegionHomologyTargetRequest {
            genome_id: format!("target_{}", role.as_str()),
            required: false,
            role,
            expected_loci: if role == gp::GenomicRegionHomologyTargetRole::ExpectedOrtholog {
                vec![gp::GenomicRegionHomologyExpectedLocus {
                    expected_locus_id: "ortholog_locus".to_string(),
                    reference: gp::GenomicRegionReference {
                        assembly_name: "target_v1".to_string(),
                        contig_name: "chrO".to_string(),
                        ..Default::default()
                    },
                    start_0based: 200,
                    end_0based_exclusive: 230,
                    evidence_id: "explicit_orthology".to_string(),
                    source_id: "synthetic_ortholog_map".to_string(),
                    ..Default::default()
                }]
            } else {
                vec![]
            },
        }
    }

    fn hsp(
        target_id: &str,
        subject: &str,
        q: &str,
        s: &str,
        q_start: usize,
        subject_start: u64,
        strand: gp::GenomicRegionStrand,
    ) -> gp::GenomicRegionHomologyHsp {
        let query_consumed = q.bytes().filter(|base| *base != b'-').count();
        let subject_consumed = s.bytes().filter(|base| *base != b'-').count() as u64;
        let (subject_min, subject_max) = if strand == gp::GenomicRegionStrand::Minus {
            (
                subject_start.saturating_sub(subject_consumed),
                subject_start,
            )
        } else {
            (subject_start, subject_start + subject_consumed)
        };
        gp::GenomicRegionHomologyHsp {
            hsp_id: short_sha256_id("hsp", &format!("{target_id}:{subject}:{q_start}:{q}:{s}")),
            target_genome_id: target_id.to_string(),
            subject_id_raw: subject.to_string(),
            subject_id: subject.to_string(),
            strand,
            query_start_0based: q_start,
            query_end_0based_exclusive: q_start + query_consumed,
            subject_start_0based: subject_min,
            subject_end_0based_exclusive: subject_max,
            identity_percent: 90.0,
            alignment_length_bp: q.len(),
            evalue: 1.0e-10,
            bit_score: 80.0,
            aligned_query: q.to_string(),
            aligned_subject: s.to_string(),
            ..Default::default()
        }
    }

    #[test]
    fn query_projection_omits_target_insertions_and_retains_them_as_provenance() {
        let query = query();
        let ortholog = target(gp::GenomicRegionHomologyTargetRole::ExpectedOrtholog);
        let effective = gp::GenomicRegionHomologyEffectiveRequest {
            set_id: query.set_id.clone(),
            region_id: query.region.region_id.clone(),
            region_content_sha256: query.region.content_sha256.clone(),
            query_genome_id: "query_genome".to_string(),
            targets: vec![ortholog.clone()],
            policy: gp::GenomicRegionHomologySearchPolicy {
                min_conserved_block_bp: 4,
                ..Default::default()
            },
            ..Default::default()
        };
        let insertion_hsp = hsp(
            &ortholog.genome_id,
            "chrO",
            "AACCGG--TTAACCGGTT",
            "AACCGGAATTAACCGGTT",
            0,
            205,
            gp::GenomicRegionStrand::Plus,
        );
        let report = finalize_projection(
            query.clone(),
            effective,
            vec![gp::GenomicRegionHomologyTargetResult {
                target: ortholog,
                status: gp::GenomicRegionHomologyTargetStatus::Available,
                accepted_hsp_count: 1,
                ..Default::default()
            }],
            vec![insertion_hsp],
            "sha256:request".to_string(),
        )
        .expect("projection");
        assert_eq!(
            report.alignment_rows[0].locus_class,
            gp::GenomicRegionHomologyLocusClass::Query
        );
        assert!(
            report
                .alignment_rows
                .iter()
                .all(|row| row.query_projection.len() == query.sequence.len())
        );
        assert_eq!(report.omitted_insertions.len(), 1);
        assert_eq!(report.omitted_insertions[0].target_sequence, "AA");
        assert_eq!(report.omitted_insertions[0].query_anchor_0based, 6);
    }

    #[test]
    fn blast_similarity_is_not_called_ortholog_without_expected_locus_overlap() {
        let query = query();
        let expected = target(gp::GenomicRegionHomologyTargetRole::ExpectedOrtholog);
        let effective = gp::GenomicRegionHomologyEffectiveRequest {
            set_id: query.set_id.clone(),
            region_id: query.region.region_id.clone(),
            region_content_sha256: query.region.content_sha256.clone(),
            query_genome_id: "query_genome".to_string(),
            targets: vec![expected.clone()],
            policy: gp::GenomicRegionHomologySearchPolicy {
                min_conserved_block_bp: 4,
                ..Default::default()
            },
            ..Default::default()
        };
        let report = finalize_projection(
            query,
            effective,
            vec![gp::GenomicRegionHomologyTargetResult {
                target: expected.clone(),
                status: gp::GenomicRegionHomologyTargetStatus::Available,
                ..Default::default()
            }],
            vec![hsp(
                &expected.genome_id,
                "chrUnexpected",
                "AACCGGTTAACCGGTT",
                "AACCGGTTAACCGGTT",
                0,
                900,
                gp::GenomicRegionStrand::Plus,
            )],
            "sha256:request".to_string(),
        )
        .expect("projection");
        assert_eq!(
            report.loci[0].locus_class,
            gp::GenomicRegionHomologyLocusClass::CrossSpeciesUnassigned
        );
        assert!(report.loci[0].orthology_evidence_id.is_none());
    }

    fn with_assessed_same_genome(
        mut report: gp::GenomicRegionHomologyScreenReport,
    ) -> gp::GenomicRegionHomologyScreenReport {
        let same = gp::GenomicRegionHomologyTargetRequest {
            genome_id: "query_genome".to_string(),
            role: gp::GenomicRegionHomologyTargetRole::SameGenome,
            ..Default::default()
        };
        report.effective_request.targets.push(same.clone());
        report.targets.push(gp::GenomicRegionHomologyTargetResult {
            target: same,
            status: gp::GenomicRegionHomologyTargetStatus::NoAcceptedSimilarity,
            ..Default::default()
        });
        resign_test_report(&mut report);
        report
    }

    fn resign_test_report(report: &mut gp::GenomicRegionHomologyScreenReport) {
        let mut content = report.clone();
        content.content_sha256.clear();
        content.op_id = None;
        content.run_id = None;
        report.content_sha256 = canonical_digest(&content, "synthetic report").expect("digest");
    }

    #[test]
    fn module_assessment_emits_traceable_standalone_and_repetitive_states() {
        let query = query();
        let expected = target(gp::GenomicRegionHomologyTargetRole::ExpectedOrtholog);
        let effective = gp::GenomicRegionHomologyEffectiveRequest {
            set_id: query.set_id.clone(),
            region_id: query.region.region_id.clone(),
            region_content_sha256: query.region.content_sha256.clone(),
            query_genome_id: "query_genome".to_string(),
            targets: vec![expected.clone()],
            policy: gp::GenomicRegionHomologySearchPolicy {
                min_conserved_block_bp: 4,
                ..Default::default()
            },
            ..Default::default()
        };
        let report = finalize_projection(
            query,
            effective,
            vec![gp::GenomicRegionHomologyTargetResult {
                target: expected.clone(),
                status: gp::GenomicRegionHomologyTargetStatus::Available,
                ..Default::default()
            }],
            vec![hsp(
                &expected.genome_id,
                "chrO",
                "AACCGGTTAACCGGTT",
                "AACCGGTTAACCGGTT",
                0,
                205,
                gp::GenomicRegionStrand::Plus,
            )],
            "sha256:request".to_string(),
        )
        .expect("projection");
        let report = with_assessed_same_genome(report);
        let engine = GentleEngine::default();
        let request = gp::PromoterModuleAssessmentRequest {
            homology_report: Box::new(report.clone()),
            selected_evidence_spans: vec![gp::PromoterModuleEvidenceSpan {
                evidence_id: "motif_tuple".to_string(),
                evidence_kind: "motif_tuple".to_string(),
                query_start_0based: 4,
                query_end_0based_exclusive: 8,
                required: true,
                source_id: "synthetic".to_string(),
                evidence_statement: "synthetic motif tuple".to_string(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let assessed = engine
            .assess_promoter_conserved_modules(request.clone(), "op", "run")
            .expect("assessment");
        assert_eq!(
            assessed.hypothesis,
            gp::PromoterModuleHypothesisKind::StandaloneReporterCandidate
        );
        assert!(!assessed.decision_trace.is_empty());

        for status in [
            None,
            Some(gp::GenomicRegionHomologyTargetStatus::Unavailable),
            Some(gp::GenomicRegionHomologyTargetStatus::SearchOutputTooBroad),
        ] {
            let mut missing = report.clone();
            if let Some(status) = status {
                missing
                    .targets
                    .last_mut()
                    .expect("same-genome result")
                    .status = status;
            } else {
                missing.targets.pop();
                missing.effective_request.targets.pop();
            }
            resign_test_report(&mut missing);
            let missing = engine
                .assess_promoter_conserved_modules(
                    gp::PromoterModuleAssessmentRequest {
                        homology_report: Box::new(missing),
                        ..request.clone()
                    },
                    "op",
                    "run",
                )
                .expect("unassessed report");
            assert_eq!(
                missing.hypothesis,
                gp::PromoterModuleHypothesisKind::InsufficientEvidence
            );
            let rule = missing
                .decision_trace
                .iter()
                .find(|rule| rule.rule_id == "same_genome_interpretation_unique")
                .expect("uniqueness rule");
            assert!(!rule.satisfied);
            assert!(rule.detail.contains("unassessed"));
        }

        let mut repetitive_report = report;
        repetitive_report.same_genome_nonself_query_coverage_percent = 100.0;
        let mut content = repetitive_report.clone();
        content.content_sha256.clear();
        content.op_id = None;
        content.run_id = None;
        repetitive_report.content_sha256 =
            canonical_digest(&content, "repetitive report").expect("digest");
        let repetitive = engine
            .assess_promoter_conserved_modules(
                gp::PromoterModuleAssessmentRequest {
                    homology_report: Box::new(repetitive_report),
                    ..request
                },
                "op",
                "run",
            )
            .expect("assessment");
        assert_eq!(
            repetitive.hypothesis,
            gp::PromoterModuleHypothesisKind::RepetitiveOrAmbiguous
        );
    }

    #[test]
    fn projection_orders_classes_and_handles_reverse_substitutions_and_deletions() {
        let query = query();
        let expected = target(gp::GenomicRegionHomologyTargetRole::ExpectedOrtholog);
        let unassigned = target(gp::GenomicRegionHomologyTargetRole::CrossSpeciesUnassigned);
        let same = gp::GenomicRegionHomologyTargetRequest {
            genome_id: "query_genome".to_string(),
            role: gp::GenomicRegionHomologyTargetRole::SameGenome,
            ..Default::default()
        };
        let effective = gp::GenomicRegionHomologyEffectiveRequest {
            set_id: query.set_id.clone(),
            region_id: query.region.region_id.clone(),
            region_content_sha256: query.region.content_sha256.clone(),
            query_genome_id: same.genome_id.clone(),
            targets: vec![same.clone(), unassigned.clone(), expected.clone()],
            policy: gp::GenomicRegionHomologySearchPolicy {
                min_conserved_block_bp: 2,
                ..Default::default()
            },
            ..Default::default()
        };
        let report = finalize_projection(
            query,
            effective,
            vec![
                gp::GenomicRegionHomologyTargetResult {
                    target: same.clone(),
                    status: gp::GenomicRegionHomologyTargetStatus::Available,
                    ..Default::default()
                },
                gp::GenomicRegionHomologyTargetResult {
                    target: unassigned.clone(),
                    status: gp::GenomicRegionHomologyTargetStatus::Available,
                    ..Default::default()
                },
                gp::GenomicRegionHomologyTargetResult {
                    target: expected.clone(),
                    status: gp::GenomicRegionHomologyTargetStatus::Available,
                    ..Default::default()
                },
            ],
            vec![
                hsp(
                    &same.genome_id,
                    "chrQ",
                    "AACCGGTTAACCGGTT",
                    "AACCGGTTAACCGGTT",
                    0,
                    100,
                    gp::GenomicRegionStrand::Plus,
                ),
                hsp(
                    &same.genome_id,
                    "chrQ",
                    "AACCGGTT",
                    "AACCGGTT",
                    0,
                    500,
                    gp::GenomicRegionStrand::Plus,
                ),
                hsp(
                    &unassigned.genome_id,
                    "chrU",
                    "AACCGGTT",
                    "AACCGGTT",
                    0,
                    300,
                    gp::GenomicRegionStrand::Plus,
                ),
                hsp(
                    &expected.genome_id,
                    "chrO",
                    "AACCGGTT",
                    "AAC-GATT",
                    0,
                    220,
                    gp::GenomicRegionStrand::Minus,
                ),
            ],
            "sha256:request".to_string(),
        )
        .expect("projection");

        assert_eq!(
            report.alignment_rows[0].locus_class,
            gp::GenomicRegionHomologyLocusClass::Query
        );
        let expected_row = report
            .alignment_rows
            .iter()
            .find(|row| row.locus_class == gp::GenomicRegionHomologyLocusClass::ExpectedOrtholog)
            .expect("expected ortholog row");
        assert_eq!(expected_row.strand, gp::GenomicRegionStrand::Minus);
        assert_eq!(&expected_row.query_projection[..8], "...-.A..");
        let classes = report
            .alignment_rows
            .iter()
            .map(|row| locus_class_rank(row.locus_class))
            .collect::<Vec<_>>();
        assert!(classes.windows(2).all(|pair| pair[0] <= pair[1]));
        assert_eq!(report.same_genome_nonself_locus_count, 1);
        assert!(report.conserved_blocks.iter().any(|block| {
            block.support_class == gp::GenomicRegionHomologySupportClass::SameGenomeNonself
        }));
    }

    #[test]
    fn hsp_chaining_and_overlap_precedence_are_deterministic() {
        let query = query();
        let unassigned = target(gp::GenomicRegionHomologyTargetRole::CrossSpeciesUnassigned);
        let effective = gp::GenomicRegionHomologyEffectiveRequest {
            set_id: query.set_id.clone(),
            region_id: query.region.region_id.clone(),
            region_content_sha256: query.region.content_sha256.clone(),
            query_genome_id: "query_genome".to_string(),
            targets: vec![unassigned.clone()],
            policy: gp::GenomicRegionHomologySearchPolicy {
                max_chain_gap_bp: 10,
                min_conserved_block_bp: 2,
                ..Default::default()
            },
            ..Default::default()
        };
        let first = hsp(
            &unassigned.genome_id,
            "chrU",
            "AACC",
            "AACC",
            0,
            100,
            gp::GenomicRegionStrand::Plus,
        );
        let second = hsp(
            &unassigned.genome_id,
            "chrU",
            "TTAA",
            "TTAA",
            6,
            106,
            gp::GenomicRegionStrand::Plus,
        );
        let mut lower = hsp(
            &unassigned.genome_id,
            "chrU",
            "AACC",
            "TTTT",
            0,
            100,
            gp::GenomicRegionStrand::Plus,
        );
        lower.hsp_id = "lower_score".to_string();
        lower.bit_score = 20.0;
        let report = finalize_projection(
            query,
            effective,
            vec![gp::GenomicRegionHomologyTargetResult {
                target: unassigned,
                status: gp::GenomicRegionHomologyTargetStatus::Available,
                ..Default::default()
            }],
            vec![second, lower, first],
            "sha256:request".to_string(),
        )
        .expect("projection");
        assert_eq!(report.loci.len(), 1);
        assert_eq!(report.loci[0].source_hsp_ids.len(), 3);
        let row = &report.alignment_rows[1];
        assert!(row.query_projection.starts_with("...."));
        assert!(!row.conflicts.is_empty());
        let serialized_once = serde_json::to_vec(&report).expect("serialize report");
        let serialized_twice = serde_json::to_vec(&report).expect("serialize report again");
        assert_eq!(serialized_once, serialized_twice);
    }

    #[test]
    fn repeated_query_segment_at_a_distant_subject_locus_is_not_chained_to_self() {
        let query = query();
        let same = target(gp::GenomicRegionHomologyTargetRole::SameGenome);
        let self_hsp = hsp(
            &same.genome_id,
            "chrQ",
            "AACC",
            "AACC",
            0,
            100,
            gp::GenomicRegionStrand::Plus,
        );
        let duplicate_hsp = hsp(
            &same.genome_id,
            "chrQ",
            "AACC",
            "AACC",
            0,
            300,
            gp::GenomicRegionStrand::Plus,
        );
        let chains = chain_hsps(&[self_hsp, duplicate_hsp], 500);
        assert_eq!(chains.len(), 2);

        let effective = gp::GenomicRegionHomologyEffectiveRequest {
            set_id: query.set_id.clone(),
            region_id: query.region.region_id.clone(),
            region_content_sha256: query.region.content_sha256.clone(),
            query_genome_id: same.genome_id.clone(),
            targets: vec![same.clone()],
            policy: gp::GenomicRegionHomologySearchPolicy {
                max_chain_gap_bp: 500,
                min_conserved_block_bp: 2,
                ..Default::default()
            },
            ..Default::default()
        };
        let report = finalize_projection(
            query,
            effective,
            vec![gp::GenomicRegionHomologyTargetResult {
                target: same,
                status: gp::GenomicRegionHomologyTargetStatus::Available,
                ..Default::default()
            }],
            chains.into_iter().flatten().collect(),
            "sha256:request".to_string(),
        )
        .expect("projection");
        assert_eq!(report.same_genome_nonself_locus_count, 1);
    }

    #[test]
    fn module_assessment_distinguishes_paired_and_insufficient_evidence() {
        let query = query();
        let expected = target(gp::GenomicRegionHomologyTargetRole::ExpectedOrtholog);
        let effective = gp::GenomicRegionHomologyEffectiveRequest {
            set_id: query.set_id.clone(),
            region_id: query.region.region_id.clone(),
            region_content_sha256: query.region.content_sha256.clone(),
            query_genome_id: "query_genome".to_string(),
            targets: vec![expected.clone()],
            policy: gp::GenomicRegionHomologySearchPolicy {
                min_conserved_block_bp: 4,
                ..Default::default()
            },
            ..Default::default()
        };
        let report = finalize_projection(
            query,
            effective,
            vec![gp::GenomicRegionHomologyTargetResult {
                target: expected.clone(),
                status: gp::GenomicRegionHomologyTargetStatus::Available,
                ..Default::default()
            }],
            vec![
                hsp(
                    &expected.genome_id,
                    "chrO",
                    "AACC",
                    "AACC",
                    0,
                    205,
                    gp::GenomicRegionStrand::Plus,
                ),
                hsp(
                    &expected.genome_id,
                    "chrO",
                    "AACC",
                    "AACC",
                    8,
                    213,
                    gp::GenomicRegionStrand::Plus,
                ),
            ],
            "sha256:request".to_string(),
        )
        .expect("projection");
        let report = with_assessed_same_genome(report);
        let engine = GentleEngine::default();
        let spans = vec![
            gp::PromoterModuleEvidenceSpan {
                evidence_id: "left".to_string(),
                query_start_0based: 0,
                query_end_0based_exclusive: 4,
                required: true,
                ..Default::default()
            },
            gp::PromoterModuleEvidenceSpan {
                evidence_id: "right".to_string(),
                query_start_0based: 8,
                query_end_0based_exclusive: 12,
                required: true,
                ..Default::default()
            },
        ];
        let paired = engine
            .assess_promoter_conserved_modules(
                gp::PromoterModuleAssessmentRequest {
                    homology_report: Box::new(report.clone()),
                    selected_evidence_spans: spans,
                    max_partner_gap_bp: 8,
                    ..Default::default()
                },
                "op",
                "run",
            )
            .expect("paired assessment");
        assert_eq!(
            paired.hypothesis,
            gp::PromoterModuleHypothesisKind::PairedContextCandidate
        );

        let insufficient = engine
            .assess_promoter_conserved_modules(
                gp::PromoterModuleAssessmentRequest {
                    homology_report: Box::new(report),
                    selected_evidence_spans: vec![gp::PromoterModuleEvidenceSpan {
                        evidence_id: "gap".to_string(),
                        query_start_0based: 5,
                        query_end_0based_exclusive: 7,
                        required: true,
                        ..Default::default()
                    }],
                    ..Default::default()
                },
                "op",
                "run",
            )
            .expect("insufficient assessment");
        assert_eq!(
            insufficient.hypothesis,
            gp::PromoterModuleHypothesisKind::InsufficientEvidence
        );
    }

    #[test]
    fn no_indexes_is_an_explicit_query_only_report_and_fingerprint_changes_cache_key() {
        let query = query();
        let effective = gp::GenomicRegionHomologyEffectiveRequest {
            set_id: query.set_id.clone(),
            region_id: query.region.region_id.clone(),
            region_content_sha256: query.region.content_sha256.clone(),
            query_genome_id: "query_genome".to_string(),
            ..Default::default()
        };
        let report = finalize_projection(
            query.clone(),
            effective.clone(),
            vec![],
            vec![],
            "sha256:request".to_string(),
        )
        .expect("query-only projection");
        assert!(report.targets.is_empty());
        assert_eq!(report.alignment_rows.len(), 1);
        assert!(
            report
                .warnings
                .iter()
                .any(|warning| warning.contains("No validated"))
        );

        let first = region_homology_cache_key(
            &query.sequence_sha256,
            &effective,
            "genome:sha256:first:blastn-1:valid:genomic_dna",
        )
        .expect("first cache key");
        let replacement = region_homology_cache_key(
            &query.sequence_sha256,
            &effective,
            "genome:sha256:replacement:blastn-1:valid:genomic_dna",
        )
        .expect("replacement cache key");
        assert_ne!(first, replacement);
    }

    #[test]
    fn subject_identifier_normalization_uses_prepared_fasta_aliases() {
        let aliases = BTreeMap::from([
            ("KI270750.1".to_string(), "KI270750.1".to_string()),
            ("ki270750.1".to_string(), "KI270750.1".to_string()),
        ]);
        assert_eq!(
            normalize_subject_id("gb|KI270750.1|", &aliases),
            "KI270750.1"
        );
    }

    #[test]
    fn optional_missing_index_is_reported_without_mutating_state_and_required_fails_preflight() {
        let temp = tempfile::tempdir().expect("temporary catalog directory");
        let catalog_path = temp.path().join("empty_genomes.json");
        std::fs::write(&catalog_path, "{}\n").expect("write empty catalog");
        let mut state = ProjectState::default();
        state.sequences.insert(
            "query_seq".to_string(),
            DNAsequence::from_sequence("AACCGGTTAACCGGTT").expect("query DNA"),
        );
        state.metadata.insert(
            PROVENANCE_METADATA_KEY.to_string(),
            serde_json::json!({
                GENOME_EXTRACTIONS_METADATA_KEY: [{
                    "seq_id": "query_seq",
                    "genome_id": "query_assembly",
                    "chromosome": "chrQ",
                    "start_1based": 101,
                    "end_1based": 116,
                    "anchor_strand": "+",
                    "anchor_verified": true,
                    "recorded_at_unix_ms": 123
                }]
            }),
        );
        let mut engine = GentleEngine::from_state(state);
        let sequence_sha256 = sha256_prefixed_bytes(
            engine
                .state()
                .sequences
                .get("query_seq")
                .expect("query sequence")
                .forward_bytes(),
        );
        engine
            .apply(Operation::CreateGenomicRegion {
                request: gp::GenomicRegionCreateRequest {
                    set_id: "set".to_string(),
                    region_id: Some("region".to_string()),
                    purpose: gp::GenomicRegionPurpose::ReporterCandidate,
                    interval: gp::GenomicRegionInterval {
                        reference: gp::GenomicRegionReference {
                            assembly_name: "query_assembly".to_string(),
                            contig_name: "chrQ".to_string(),
                            ..Default::default()
                        },
                        start_0based: 100,
                        end_0based_exclusive: 116,
                        strand: gp::GenomicRegionStrand::Plus,
                        ..Default::default()
                    },
                    local_projection: Some(gp::GenomicRegionLocalProjection {
                        seq_id: "query_seq".to_string(),
                        sequence_sha256,
                        source_genome_id: "query_assembly".to_string(),
                        anchor_start_1based: 101,
                        anchor_end_1based: 116,
                        anchor_strand: gp::GenomicRegionStrand::Plus,
                        local_start_0based: 0,
                        local_end_0based_exclusive: 16,
                        local_strand: gp::GenomicRegionStrand::Plus,
                        status: gp::GenomicRegionLocalProjectionStatus::Current,
                    }),
                    ..Default::default()
                },
            })
            .expect("capture query region");
        let state_before = serde_json::to_vec(engine.state()).expect("serialize state");
        let base_request = gp::GenomicRegionHomologyScreenRequest {
            set_id: "set".to_string(),
            region_id: "region".to_string(),
            query_genome_id: Some("query_assembly".to_string()),
            catalog_path: Some(catalog_path.display().to_string()),
            targets: vec![gp::GenomicRegionHomologyTargetRequest {
                genome_id: "missing_target".to_string(),
                required: false,
                ..Default::default()
            }],
            ..Default::default()
        };
        let report = engine
            .screen_genomic_region_homology(base_request.clone(), "op", "run", &mut |_| true)
            .expect("optional missing target remains reportable");
        assert_eq!(
            report.targets[0].status,
            gp::GenomicRegionHomologyTargetStatus::Unavailable
        );
        assert_eq!(
            state_before,
            serde_json::to_vec(engine.state()).expect("serialize state after screen")
        );

        let error = engine
            .screen_genomic_region_homology(
                gp::GenomicRegionHomologyScreenRequest {
                    targets: vec![gp::GenomicRegionHomologyTargetRequest {
                        required: true,
                        ..base_request.targets[0].clone()
                    }],
                    ..base_request
                },
                "op",
                "run",
                &mut |_| true,
            )
            .expect_err("required missing target must fail before BLAST");
        assert!(error.message.contains("required homology target"));
    }

    #[test]
    fn opt_in_real_blast_region_homology_smoke() {
        if std::env::var_os("GENTLE_TEST_REGION_HOMOLOGY_BLAST").is_none() {
            return;
        }
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let asset_dir = root.join("docs/examples/assets/region_homology_demo");
        let catalog_path = asset_dir.join("genomes.json");
        let cache = tempfile::tempdir().expect("temporary prepared-genome cache");
        let cache_path = cache.path().display().to_string();
        let catalog_string = catalog_path.display().to_string();
        let mut engine = GentleEngine::from_state(ProjectState::default());
        for genome_id in [
            "QueryRegionToy",
            "ExpectedOrthologToy",
            "UnassignedSimilarityToy",
        ] {
            engine
                .apply(Operation::PrepareGenome {
                    genome_id: genome_id.to_string(),
                    catalog_path: Some(catalog_string.clone()),
                    cache_dir: Some(cache_path.clone()),
                    timeout_seconds: Some(60),
                })
                .unwrap_or_else(|error| panic!("prepare {genome_id}: {}", error.message));
        }
        engine
            .apply(Operation::ExtractGenomeRegion {
                genome_id: "QueryRegionToy".to_string(),
                chromosome: "chrQ".to_string(),
                start_1based: 1,
                end_1based: 188,
                output_id: Some("query_seq".to_string()),
                annotation_scope: Some(GenomeAnnotationScope::None),
                max_annotation_features: None,
                include_genomic_annotation: Some(false),
                catalog_path: Some(catalog_string.clone()),
                cache_dir: Some(cache_path.clone()),
            })
            .expect("extract anchored synthetic query");
        engine
            .apply(Operation::CaptureGenomicRegion {
                request: gp::GenomicRegionCaptureRequest {
                    set_id: "set".to_string(),
                    region_id: Some("query".to_string()),
                    purpose: gp::GenomicRegionPurpose::ReporterCandidate,
                    source: gp::GenomicRegionCaptureSource::SequenceSelection {
                        seq_id: "query_seq".to_string(),
                        local_start_0based: 20,
                        local_end_0based_exclusive: 100,
                        strand: gp::GenomicRegionStrand::Plus,
                        reference_override: Some(gp::GenomicRegionReference {
                            assembly_name: "QueryRegionToy".to_string(),
                            contig_name: "chrQ".to_string(),
                            ..Default::default()
                        }),
                    },
                    ..Default::default()
                },
            })
            .expect("capture synthetic query");
        let report = engine
            .screen_genomic_region_homology(
                gp::GenomicRegionHomologyScreenRequest {
                    set_id: "set".to_string(),
                    region_id: "query".to_string(),
                    query_genome_id: Some("QueryRegionToy".to_string()),
                    catalog_path: Some(catalog_string),
                    cache_dir: Some(cache_path),
                    targets: vec![
                        gp::GenomicRegionHomologyTargetRequest {
                            genome_id: "QueryRegionToy".to_string(),
                            required: true,
                            role: gp::GenomicRegionHomologyTargetRole::SameGenome,
                            ..Default::default()
                        },
                        gp::GenomicRegionHomologyTargetRequest {
                            genome_id: "ExpectedOrthologToy".to_string(),
                            required: true,
                            role: gp::GenomicRegionHomologyTargetRole::ExpectedOrtholog,
                            expected_loci: vec![gp::GenomicRegionHomologyExpectedLocus {
                                expected_locus_id: "expected".to_string(),
                                reference: gp::GenomicRegionReference {
                                    assembly_name: "ExpectedOrthologToy".to_string(),
                                    contig_name: "chrO".to_string(),
                                    ..Default::default()
                                },
                                start_0based: 30,
                                end_0based_exclusive: 111,
                                strand: gp::GenomicRegionStrand::Plus,
                                evidence_id: "synthetic_expected_locus".to_string(),
                                source_id: "region_homology_demo_v1".to_string(),
                                ..Default::default()
                            }],
                        },
                        gp::GenomicRegionHomologyTargetRequest {
                            genome_id: "UnassignedSimilarityToy".to_string(),
                            role: gp::GenomicRegionHomologyTargetRole::CrossSpeciesUnassigned,
                            ..Default::default()
                        },
                    ],
                    policy: gp::GenomicRegionHomologySearchPolicy {
                        min_alignment_length_bp: 12,
                        min_conserved_block_bp: 12,
                        max_hsps_per_target: 1_000,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                "op",
                "run",
                &mut |_| true,
            )
            .expect("real BLAST homology screen");
        assert!(
            report
                .alignment_rows
                .iter()
                .all(|row| { row.query_projection.len() == report.query.sequence.len() })
        );
        assert!(
            report.loci.iter().any(|locus| {
                locus.locus_class == gp::GenomicRegionHomologyLocusClass::ExpectedOrtholog
            }),
            "explicit expected locus was not recognized: {:#?}",
            report.loci
        );
        assert!(report.same_genome_nonself_locus_count > 0);
        assert!(!report.conserved_blocks.is_empty());
    }
}
