//! Deterministic, offline screening of primer/probe binding sites against VCF.
//!
//! Candidate geometry is explicit and assembly-aware. The screen never treats
//! source catalogue claims as validation, never assigns historical coordinates
//! to a current reference, and scans one local VCF once for all supplied pairs.

use crate::{
    digest_utils::{
        canonical_oligo_sequence, oligo_full_id, primer_pair_full_id, sha256_file_hex,
        sha256_prefixed_bytes,
    },
    engine::{
        PRIMER_VARIANT_EVIDENCE_SCHEMA, PRIMER_VARIANT_RESOURCE_MANIFEST_SCHEMA,
        PRIMER_VARIANT_SCREEN_REQUEST_SCHEMA, PRIMER_VARIANT_SCREEN_SCHEMA,
        PrimerVariantAmpliconSegment, PrimerVariantBindingSegment, PrimerVariantBindingStrand,
        PrimerVariantCandidateSource, PrimerVariantEvidenceReport, PrimerVariantEvidenceStatus,
        PrimerVariantKind, PrimerVariantOverlap, PrimerVariantOverlapKind,
        PrimerVariantProbeOverlapPolicy, PrimerVariantResourceManifest,
        PrimerVariantScreenCandidate, PrimerVariantScreenOligo, PrimerVariantScreenReport,
        PrimerVariantScreenRequest, PrimerVariantSourceInput,
    },
};
use flate2::read::MultiGzDecoder;
use serde::Serialize;
use std::{
    collections::{BTreeMap, BTreeSet, HashSet},
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

#[derive(Debug, Clone)]
struct ResolvedVariantSource {
    vcf_path: PathBuf,
    reference_assembly: String,
    source_name: String,
    source_release: String,
    population: String,
    retrieval_time: String,
    allele_frequency_info_field: Option<String>,
    expected_content_sha256: Option<String>,
    contig_aliases: BTreeMap<String, String>,
    manifest_sha256: Option<String>,
}

#[derive(Debug, Clone)]
struct PairWork {
    pair_id: String,
    geometry_sha256: String,
    forward: PrimerVariantScreenOligo,
    reverse: PrimerVariantScreenOligo,
    probe: Option<PrimerVariantScreenOligo>,
    amplicon_segments: Vec<PrimerVariantAmpliconSegment>,
    sources: Vec<PrimerVariantCandidateSource>,
    reference_names: BTreeSet<String>,
    overlaps: Vec<PrimerVariantOverlap>,
    overlapping_record_numbers: BTreeSet<usize>,
    reference_incompatible: bool,
    warnings: Vec<String>,
}

#[derive(Debug, Clone)]
struct ParsedVcfRecord {
    reference_name: String,
    position_1based: usize,
    id: String,
    reference_allele: String,
    alternate_alleles: Vec<String>,
    filter: String,
    info: String,
}

#[derive(Debug, Clone)]
struct OligoOverlap {
    oligo_position_1based: usize,
    oligo_end_position_1based: usize,
    distance_from_3prime_end: usize,
    critical_three_prime: bool,
    reference_match: Option<bool>,
}

/// Screen every unique physical pair in one request against one local VCF.
pub fn screen_primer_variants(
    request: PrimerVariantScreenRequest,
) -> Result<PrimerVariantScreenReport, String> {
    if request.schema != PRIMER_VARIANT_SCREEN_REQUEST_SCHEMA {
        return Err(format!(
            "Unsupported primer variant screen request schema '{}'; expected '{}'",
            request.schema, PRIMER_VARIANT_SCREEN_REQUEST_SCHEMA
        ));
    }
    let reference_assembly = required_text(&request.reference_assembly, "reference_assembly")?;
    if request.candidates.is_empty() {
        return Err("Primer variant screen requires at least one candidate".to_string());
    }
    if request.critical_3prime_bases == 0 {
        return Err("critical_3prime_bases must be >= 1".to_string());
    }
    if let Some(value) = request.maximum_allowed_frequency
        && !(0.0..=1.0).contains(&value)
    {
        return Err(format!(
            "maximum_allowed_frequency must be between 0 and 1, got {value}"
        ));
    }

    let source = resolve_variant_source(&request.source)?;
    let content_sha256 = format!(
        "sha256:{}",
        sha256_file_hex(&source.vcf_path).map_err(|error| {
            format!(
                "Could not fingerprint variant VCF '{}': {error}",
                source.vcf_path.display()
            )
        })?
    );
    if let Some(expected) = source.expected_content_sha256.as_deref() {
        let expected = normalize_sha256(expected);
        if expected != content_sha256 {
            return Err(format!(
                "Variant VCF '{}' has digest '{}', expected '{}'",
                source.vcf_path.display(),
                content_sha256,
                expected
            ));
        }
    }

    let mut pairs = normalize_candidates(&request.candidates, &source.contig_aliases)?;
    let request_sha256 = fingerprint_request(&request, &source, &content_sha256)?;
    let screen_id = short_digest_id("primer_variant_screen", &request_sha256);
    let assembly_compatible = assembly_matches(&reference_assembly, &source.reference_assembly);
    if !assembly_compatible {
        let warning = format!(
            "Candidate coordinates use assembly '{}' but the variant source declares '{}'.",
            reference_assembly, source.reference_assembly
        );
        for pair in pairs.values_mut() {
            pair.reference_incompatible = true;
            pair.warnings.push(warning.clone());
        }
    }

    let mut global_warnings = Vec::new();
    let mut vcf_record_count = 0usize;
    let mut globally_overlapping_records = BTreeSet::new();
    let mut declared_contigs = BTreeSet::new();
    let mut observed_contigs = BTreeSet::new();
    let mut declared_reference = None;
    if assembly_compatible {
        visit_vcf_records(
            &source.vcf_path,
            |header| {
                if let Some(value) = header.strip_prefix("##reference=") {
                    declared_reference = Some(value.trim().to_string());
                }
                if let Some(contig) = parse_vcf_contig_header(header) {
                    declared_contigs.insert(contig);
                }
            },
            |record_number, record| {
                vcf_record_count = vcf_record_count.saturating_add(1);
                observed_contigs.insert(record.reference_name.clone());
                if !pairs
                    .values()
                    .any(|pair| pair_may_overlap_record(pair, &record, &source.contig_aliases))
                {
                    return Ok(());
                }
                let frequencies = allele_frequencies(
                    &record.info,
                    source.allele_frequency_info_field.as_deref(),
                    record.alternate_alleles.len(),
                    &mut global_warnings,
                    &record,
                );
                let mut record_overlapped = false;
                for pair in pairs.values_mut() {
                    let before = pair.overlaps.len();
                    screen_record_for_pair(
                        pair,
                        &record,
                        &frequencies,
                        request.maximum_allowed_frequency,
                        request.critical_3prime_bases,
                        request.probe_overlap_policy,
                        &source.contig_aliases,
                    );
                    if pair.overlaps.len() > before {
                        pair.overlapping_record_numbers.insert(record_number);
                        record_overlapped = true;
                    }
                }
                if record_overlapped {
                    globally_overlapping_records.insert(record_number);
                }
                Ok(())
            },
        )?;
    }

    if let Some(declared) = declared_reference.as_deref()
        && !declared
            .to_ascii_lowercase()
            .contains(&source.reference_assembly.to_ascii_lowercase())
    {
        let warning = format!(
            "VCF header reference '{}' does not identify declared assembly '{}'.",
            declared, source.reference_assembly
        );
        for pair in pairs.values_mut() {
            pair.reference_incompatible = true;
            pair.warnings.push(warning.clone());
        }
    }

    if assembly_compatible {
        let available_contigs = if declared_contigs.is_empty() {
            &observed_contigs
        } else {
            &declared_contigs
        };
        let available_keys = available_contigs
            .iter()
            .map(|name| contig_key(name, &source.contig_aliases))
            .collect::<BTreeSet<_>>();
        for pair in pairs.values_mut() {
            let has_available_contig = pair
                .reference_names
                .iter()
                .any(|name| available_keys.contains(&contig_key(name, &source.contig_aliases)));
            if !has_available_contig {
                pair.reference_incompatible = true;
                let evidence_kind = if declared_contigs.is_empty() {
                    "observed records"
                } else {
                    "contig declarations"
                };
                pair.warnings.push(format!(
                    "VCF {evidence_kind} do not contain any target reference for pair '{}'.",
                    pair.pair_id
                ));
            }
        }
    }

    global_warnings.sort();
    global_warnings.dedup();
    let evidence_reports = pairs
        .into_values()
        .map(|mut pair| {
            pair.overlaps.sort_by(|left, right| {
                left.reference_name
                    .cmp(&right.reference_name)
                    .then(left.position_1based.cmp(&right.position_1based))
                    .then(left.role.cmp(&right.role))
                    .then(left.alternate_allele.cmp(&right.alternate_allele))
            });
            pair.sources.sort_by(|left, right| {
                left.candidate_id
                    .cmp(&right.candidate_id)
                    .then(left.source_id.cmp(&right.source_id))
            });
            pair.warnings.extend(global_warnings.iter().cloned());
            pair.warnings.sort();
            pair.warnings.dedup();
            let status = if pair.reference_incompatible
                || pair
                    .overlaps
                    .iter()
                    .any(|overlap| overlap.reference_match != Some(true))
            {
                PrimerVariantEvidenceStatus::IncompatibleReference
            } else if pair
                .overlaps
                .iter()
                .any(|overlap| overlap.relevant_under_policy)
            {
                PrimerVariantEvidenceStatus::VariantDetected
            } else {
                PrimerVariantEvidenceStatus::EvaluatedNoRelevantVariant
            };
            let evidence_identity = serde_json::to_vec(&(
                &pair.pair_id,
                &request_sha256,
                &content_sha256,
                &pair.geometry_sha256,
            ))
            .map_err(|error| format!("Could not identify primer variant evidence: {error}"))?;
            Ok(PrimerVariantEvidenceReport {
                schema: PRIMER_VARIANT_EVIDENCE_SCHEMA.to_string(),
                evidence_id: short_digest_id(
                    "primer_variant_evidence",
                    &sha256_prefixed_bytes(&evidence_identity),
                ),
                pair_id: pair.pair_id,
                reference_assembly: reference_assembly.clone(),
                source_name: source.source_name.clone(),
                source_release: source.source_release.clone(),
                population: source.population.clone(),
                maximum_allowed_frequency: request.maximum_allowed_frequency,
                normalization_method: "vcf_record_span_overlap_without_variant_renormalization_v1"
                    .to_string(),
                retrieval_time: source.retrieval_time.clone(),
                content_sha256: content_sha256.clone(),
                request_sha256: request_sha256.clone(),
                source_path: source.vcf_path.display().to_string(),
                candidate_sources: pair.sources,
                critical_3prime_bases: request.critical_3prime_bases,
                probe_overlap_policy: request.probe_overlap_policy,
                screened_reference_names: pair.reference_names.into_iter().collect(),
                vcf_record_count,
                overlapping_record_count: pair.overlapping_record_numbers.len(),
                status,
                overlaps: pair.overlaps,
                warnings: pair.warnings,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;

    Ok(PrimerVariantScreenReport {
        schema: PRIMER_VARIANT_SCREEN_SCHEMA.to_string(),
        screen_id,
        request_sha256,
        reference_assembly,
        source_name: source.source_name,
        source_release: source.source_release,
        population: source.population,
        source_path: source.vcf_path.display().to_string(),
        content_sha256,
        candidate_count: request.candidates.len(),
        unique_pair_count: evidence_reports.len(),
        vcf_record_count,
        overlapping_record_count: globally_overlapping_records.len(),
        evidence_reports,
        warnings: global_warnings,
    })
}

fn resolve_variant_source(
    input: &PrimerVariantSourceInput,
) -> Result<ResolvedVariantSource, String> {
    match (&input.vcf_path, &input.manifest_path) {
        (Some(_), Some(_)) => {
            return Err(
                "Primer variant source accepts vcf_path or manifest_path, not both".to_string(),
            );
        }
        (None, None) => {
            return Err("Primer variant source requires vcf_path or manifest_path".to_string());
        }
        (Some(path), None) => {
            let resolved = ResolvedVariantSource {
                vcf_path: PathBuf::from(required_text(path, "source.vcf_path")?),
                reference_assembly: required_text(
                    &input.reference_assembly,
                    "source.reference_assembly",
                )?,
                source_name: required_text(&input.source_name, "source.source_name")?,
                source_release: required_text(&input.source_release, "source.source_release")?,
                population: required_text(&input.population, "source.population")?,
                retrieval_time: required_text(&input.retrieval_time, "source.retrieval_time")?,
                allele_frequency_info_field: normalize_optional_text(
                    input.allele_frequency_info_field.as_deref(),
                ),
                expected_content_sha256: normalize_optional_text(
                    input.expected_content_sha256.as_deref(),
                ),
                contig_aliases: input.contig_aliases.clone(),
                manifest_sha256: None,
            };
            ensure_vcf_exists(&resolved.vcf_path)?;
            Ok(resolved)
        }
        (None, Some(manifest_path)) => {
            let manifest_path =
                PathBuf::from(required_text(manifest_path, "source.manifest_path")?);
            let bytes = std::fs::read(&manifest_path).map_err(|error| {
                format!(
                    "Could not read primer variant resource manifest '{}': {error}",
                    manifest_path.display()
                )
            })?;
            let manifest: PrimerVariantResourceManifest =
                serde_json::from_slice(&bytes).map_err(|error| {
                    format!(
                        "Could not parse primer variant resource manifest '{}': {error}",
                        manifest_path.display()
                    )
                })?;
            if manifest.schema != PRIMER_VARIANT_RESOURCE_MANIFEST_SCHEMA {
                return Err(format!(
                    "Unsupported primer variant resource manifest schema '{}'; expected '{}'",
                    manifest.schema, PRIMER_VARIANT_RESOURCE_MANIFEST_SCHEMA
                ));
            }
            require_matching_override(
                &input.reference_assembly,
                &manifest.reference_assembly,
                "source.reference_assembly",
            )?;
            require_matching_override(
                &input.source_name,
                &manifest.source_name,
                "source.source_name",
            )?;
            require_matching_override(
                &input.source_release,
                &manifest.source_release,
                "source.source_release",
            )?;
            require_matching_override(
                &input.population,
                &manifest.population,
                "source.population",
            )?;
            require_matching_override(
                &input.retrieval_time,
                &manifest.retrieval_time,
                "source.retrieval_time",
            )?;
            require_matching_optional_override(
                input.allele_frequency_info_field.as_deref(),
                manifest.allele_frequency_info_field.as_deref(),
                "source.allele_frequency_info_field",
            )?;
            require_matching_optional_override(
                input.expected_content_sha256.as_deref(),
                manifest.expected_content_sha256.as_deref(),
                "source.expected_content_sha256",
            )?;
            let mut aliases = manifest.contig_aliases.clone();
            for (alias, canonical) in &input.contig_aliases {
                if let Some(existing) = aliases.get(alias)
                    && !existing.eq_ignore_ascii_case(canonical)
                {
                    return Err(format!(
                        "Conflicting contig alias '{}' in request ('{}') and manifest ('{}')",
                        alias, canonical, existing
                    ));
                }
                aliases.insert(alias.clone(), canonical.clone());
            }
            let raw_vcf_path =
                PathBuf::from(required_text(&manifest.vcf_path, "manifest.vcf_path")?);
            let vcf_path = if raw_vcf_path.is_absolute() {
                raw_vcf_path
            } else {
                manifest_path
                    .parent()
                    .unwrap_or_else(|| Path::new("."))
                    .join(raw_vcf_path)
            };
            ensure_vcf_exists(&vcf_path)?;
            Ok(ResolvedVariantSource {
                vcf_path,
                reference_assembly: required_text(
                    &manifest.reference_assembly,
                    "manifest.reference_assembly",
                )?,
                source_name: required_text(&manifest.source_name, "manifest.source_name")?,
                source_release: required_text(&manifest.source_release, "manifest.source_release")?,
                population: required_text(&manifest.population, "manifest.population")?,
                retrieval_time: required_text(&manifest.retrieval_time, "manifest.retrieval_time")?,
                allele_frequency_info_field: normalize_optional_text(
                    manifest.allele_frequency_info_field.as_deref(),
                ),
                expected_content_sha256: normalize_optional_text(
                    manifest.expected_content_sha256.as_deref(),
                ),
                contig_aliases: aliases,
                manifest_sha256: Some(sha256_prefixed_bytes(&bytes)),
            })
        }
    }
}

fn normalize_candidates(
    candidates: &[PrimerVariantScreenCandidate],
    aliases: &BTreeMap<String, String>,
) -> Result<BTreeMap<String, PairWork>, String> {
    let mut pairs = BTreeMap::new();
    let mut candidate_ids = HashSet::new();
    for candidate in candidates {
        let candidate_id = required_text(&candidate.candidate_id, "candidate_id")?;
        if !candidate_ids.insert(candidate_id.clone()) {
            return Err(format!("Duplicate candidate_id '{candidate_id}'"));
        }
        let forward = normalize_oligo("forward", &candidate.forward)?;
        let reverse = normalize_oligo("reverse", &candidate.reverse)?;
        let probe = candidate
            .probe
            .as_ref()
            .map(|oligo| normalize_oligo("probe", oligo))
            .transpose()?;
        let amplicon_segments = candidate
            .amplicon_segments
            .iter()
            .map(normalize_amplicon_segment)
            .collect::<Result<Vec<_>, _>>()?;
        let expected_pair_id =
            primer_pair_full_id(&forward.sequence_5_to_3, &reverse.sequence_5_to_3);
        if let Some(provided) = candidate
            .pair_id
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            && provided != expected_pair_id
        {
            return Err(format!(
                "Candidate '{}' pair_id '{}' does not match canonical physical identity '{}'",
                candidate_id, provided, expected_pair_id
            ));
        }
        let reference_names = forward
            .binding_segments
            .iter()
            .chain(reverse.binding_segments.iter())
            .chain(probe.iter().flat_map(|oligo| oligo.binding_segments.iter()))
            .map(|segment| segment.reference_name.clone())
            .chain(
                amplicon_segments
                    .iter()
                    .map(|segment| segment.reference_name.clone()),
            )
            .collect::<BTreeSet<_>>();
        let reference_keys = reference_names
            .iter()
            .map(|name| contig_key(name, aliases))
            .collect::<BTreeSet<_>>();
        if reference_keys.len() != 1 {
            return Err(format!(
                "Candidate '{}' spans multiple reference contigs: {:?}",
                candidate_id, reference_names
            ));
        }
        let geometry_sha256 =
            fingerprint_geometry(&forward, &reverse, probe.as_ref(), &amplicon_segments)?;
        let mut source = candidate.source.clone();
        source.candidate_id = candidate_id;
        match pairs.entry(expected_pair_id.clone()) {
            std::collections::btree_map::Entry::Vacant(entry) => {
                entry.insert(PairWork {
                    pair_id: expected_pair_id,
                    geometry_sha256,
                    forward,
                    reverse,
                    probe,
                    amplicon_segments,
                    sources: vec![source],
                    reference_names,
                    overlaps: Vec::new(),
                    overlapping_record_numbers: BTreeSet::new(),
                    reference_incompatible: false,
                    warnings: Vec::new(),
                });
            }
            std::collections::btree_map::Entry::Occupied(mut entry) => {
                if entry.get().geometry_sha256 != geometry_sha256 {
                    return Err(format!(
                        "Candidates sharing physical pair '{}' declare different binding geometry",
                        expected_pair_id
                    ));
                }
                entry.get_mut().sources.push(source);
            }
        }
    }
    Ok(pairs)
}

fn normalize_oligo(
    role: &str,
    oligo: &PrimerVariantScreenOligo,
) -> Result<PrimerVariantScreenOligo, String> {
    let sequence = normalize_dna(&oligo.sequence_5_to_3, &format!("{role} sequence"))?;
    if oligo.binding_segments.is_empty() {
        return Err(format!(
            "{role} primer requires at least one binding segment"
        ));
    }
    let mut segments = oligo
        .binding_segments
        .iter()
        .map(|segment| normalize_binding_segment(role, segment, &sequence))
        .collect::<Result<Vec<_>, _>>()?;
    segments.sort_by_key(|segment| segment.oligo_start_0based);
    for window in segments.windows(2) {
        if window[0].oligo_end_0based_exclusive > window[1].oligo_start_0based {
            return Err(format!(
                "{role} binding segments overlap in oligo coordinates"
            ));
        }
    }
    if segments
        .last()
        .is_none_or(|segment| segment.oligo_end_0based_exclusive != sequence.len())
    {
        return Err(format!(
            "{role} binding geometry must include the oligo 3-prime end at position {}",
            sequence.len()
        ));
    }
    Ok(PrimerVariantScreenOligo {
        sequence_5_to_3: sequence,
        binding_segments: segments,
    })
}

fn normalize_binding_segment(
    role: &str,
    segment: &PrimerVariantBindingSegment,
    oligo_sequence: &str,
) -> Result<PrimerVariantBindingSegment, String> {
    let reference_name = required_text(&segment.reference_name, "binding reference_name")?;
    validate_closed_range(
        segment.start_1based,
        segment.end_1based,
        &format!("{role} binding segment"),
    )?;
    if segment.oligo_end_0based_exclusive <= segment.oligo_start_0based
        || segment.oligo_end_0based_exclusive > oligo_sequence.len()
    {
        return Err(format!(
            "{role} binding segment has invalid oligo range {}..{} for length {}",
            segment.oligo_start_0based,
            segment.oligo_end_0based_exclusive,
            oligo_sequence.len()
        ));
    }
    let genomic_len = segment.end_1based - segment.start_1based + 1;
    let oligo_len = segment.oligo_end_0based_exclusive - segment.oligo_start_0based;
    if genomic_len != oligo_len {
        return Err(format!(
            "{role} binding segment genomic length {genomic_len} differs from oligo span {oligo_len}"
        ));
    }
    let reference_sequence = normalize_dna(
        &segment.reference_sequence_5_to_3,
        &format!("{role} binding reference sequence"),
    )?;
    if reference_sequence.len() != genomic_len {
        return Err(format!(
            "{role} binding reference sequence length {} differs from genomic span {genomic_len}",
            reference_sequence.len()
        ));
    }
    let expected_oligo_slice =
        &oligo_sequence[segment.oligo_start_0based..segment.oligo_end_0based_exclusive];
    let aligned_reference = match segment.strand {
        PrimerVariantBindingStrand::Plus => reference_sequence.clone(),
        PrimerVariantBindingStrand::Minus => reverse_complement(&reference_sequence)?,
    };
    if aligned_reference != expected_oligo_slice {
        return Err(format!(
            "{role} sequence does not match its declared {:?}-strand reference segment {}:{}-{}",
            segment.strand, reference_name, segment.start_1based, segment.end_1based
        ));
    }
    Ok(PrimerVariantBindingSegment {
        reference_name,
        start_1based: segment.start_1based,
        end_1based: segment.end_1based,
        strand: segment.strand,
        oligo_start_0based: segment.oligo_start_0based,
        oligo_end_0based_exclusive: segment.oligo_end_0based_exclusive,
        reference_sequence_5_to_3: reference_sequence,
    })
}

fn normalize_amplicon_segment(
    segment: &PrimerVariantAmpliconSegment,
) -> Result<PrimerVariantAmpliconSegment, String> {
    let reference_name = required_text(&segment.reference_name, "amplicon reference_name")?;
    validate_closed_range(segment.start_1based, segment.end_1based, "amplicon segment")?;
    let sequence = normalize_dna(
        &segment.reference_sequence_5_to_3,
        "amplicon reference sequence",
    )?;
    let expected_len = segment.end_1based - segment.start_1based + 1;
    if sequence.len() != expected_len {
        return Err(format!(
            "Amplicon reference sequence length {} differs from genomic span {expected_len}",
            sequence.len()
        ));
    }
    Ok(PrimerVariantAmpliconSegment {
        reference_name,
        start_1based: segment.start_1based,
        end_1based: segment.end_1based,
        reference_sequence_5_to_3: sequence,
    })
}

fn screen_record_for_pair(
    pair: &mut PairWork,
    record: &ParsedVcfRecord,
    frequencies: &[Option<f64>],
    maximum_allowed_frequency: Option<f64>,
    critical_3prime_bases: usize,
    probe_policy: PrimerVariantProbeOverlapPolicy,
    aliases: &BTreeMap<String, String>,
) {
    let record_end = record
        .position_1based
        .saturating_add(record.reference_allele.len().saturating_sub(1));
    for (alternate_index, alternate) in record.alternate_alleles.iter().enumerate() {
        let allele_frequency = frequencies.get(alternate_index).copied().flatten();
        let variant_id = variant_id(record, alternate, alternate_index);
        let variant_kind = classify_variant(&record.reference_allele, alternate);
        let mut oligo_rows = Vec::new();
        if let Some(hit) = oligo_overlap(&pair.forward, record, critical_3prime_bases, aliases) {
            oligo_rows.push(("forward", PrimerVariantOverlapKind::Primer, hit));
        }
        if let Some(hit) = oligo_overlap(&pair.reverse, record, critical_3prime_bases, aliases) {
            oligo_rows.push(("reverse", PrimerVariantOverlapKind::Primer, hit));
        }
        let probe_hit = pair
            .probe
            .as_ref()
            .and_then(|probe| oligo_overlap(probe, record, critical_3prime_bases, aliases));
        if probe_policy != PrimerVariantProbeOverlapPolicy::Ignore
            && let Some(hit) = probe_hit.clone()
        {
            oligo_rows.push(("probe", PrimerVariantOverlapKind::Probe, hit));
        }
        if !oligo_rows.is_empty() {
            for (role, overlap_kind, hit) in oligo_rows {
                let frequency_relevant = maximum_allowed_frequency
                    .map(|threshold| allele_frequency.is_none_or(|value| value > threshold))
                    .unwrap_or(true);
                let relevant_under_policy = frequency_relevant
                    && (overlap_kind == PrimerVariantOverlapKind::Primer
                        || probe_policy == PrimerVariantProbeOverlapPolicy::Relevant);
                if hit.reference_match != Some(true) {
                    pair.reference_incompatible = true;
                }
                if matches!(
                    variant_kind,
                    PrimerVariantKind::Insertion
                        | PrimerVariantKind::Deletion
                        | PrimerVariantKind::Complex
                ) {
                    pair.warnings.push(format!(
                        "Variant '{}' is {:?}; overlap is reported conservatively without haplotype realignment.",
                        variant_id, variant_kind
                    ));
                }
                pair.overlaps.push(PrimerVariantOverlap {
                    oligo_id: match role {
                        "forward" => oligo_full_id(&pair.forward.sequence_5_to_3),
                        "reverse" => oligo_full_id(&pair.reverse.sequence_5_to_3),
                        _ => pair
                            .probe
                            .as_ref()
                            .map(|probe| oligo_full_id(&probe.sequence_5_to_3))
                            .unwrap_or_default(),
                    },
                    role: role.to_string(),
                    variant_id: variant_id.clone(),
                    reference_name: record.reference_name.clone(),
                    position_1based: record.position_1based,
                    reference_end_1based: record_end,
                    reference_allele: record.reference_allele.clone(),
                    alternate_allele: alternate.clone(),
                    overlap_kind,
                    variant_kind,
                    oligo_position_1based: Some(hit.oligo_position_1based),
                    oligo_end_position_1based: Some(hit.oligo_end_position_1based),
                    distance_from_3prime_end: Some(hit.distance_from_3prime_end),
                    critical_three_prime: hit.critical_three_prime,
                    reference_match: hit.reference_match,
                    source_filter: record.filter.clone(),
                    allele_frequency,
                    relevant_under_policy,
                    note: overlap_note(
                        overlap_kind,
                        variant_kind,
                        allele_frequency,
                        maximum_allowed_frequency,
                        hit.critical_three_prime,
                        hit.reference_match,
                    ),
                });
            }
            continue;
        }

        if probe_policy == PrimerVariantProbeOverlapPolicy::Ignore && probe_hit.is_some() {
            continue;
        }

        if let Some(reference_match) = amplicon_overlap(pair, record, aliases) {
            if reference_match != Some(true) {
                pair.reference_incompatible = true;
            }
            pair.overlaps.push(PrimerVariantOverlap {
                oligo_id: String::new(),
                role: "amplicon".to_string(),
                variant_id,
                reference_name: record.reference_name.clone(),
                position_1based: record.position_1based,
                reference_end_1based: record_end,
                reference_allele: record.reference_allele.clone(),
                alternate_allele: alternate.clone(),
                overlap_kind: PrimerVariantOverlapKind::AmpliconOnly,
                variant_kind,
                reference_match,
                source_filter: record.filter.clone(),
                allele_frequency,
                relevant_under_policy: false,
                note: overlap_note(
                    PrimerVariantOverlapKind::AmpliconOnly,
                    variant_kind,
                    allele_frequency,
                    maximum_allowed_frequency,
                    false,
                    reference_match,
                ),
                ..PrimerVariantOverlap::default()
            });
        }
    }
}

fn pair_may_overlap_record(
    pair: &PairWork,
    record: &ParsedVcfRecord,
    aliases: &BTreeMap<String, String>,
) -> bool {
    let record_start = record.position_1based;
    let record_end = record_start.saturating_add(record.reference_allele.len().saturating_sub(1));
    let record_key = contig_key(&record.reference_name, aliases);
    pair.forward
        .binding_segments
        .iter()
        .chain(pair.reverse.binding_segments.iter())
        .chain(
            pair.probe
                .iter()
                .flat_map(|probe| probe.binding_segments.iter()),
        )
        .any(|segment| {
            contig_key(&segment.reference_name, aliases) == record_key
                && record_end >= segment.start_1based
                && record_start <= segment.end_1based
        })
        || pair.amplicon_segments.iter().any(|segment| {
            contig_key(&segment.reference_name, aliases) == record_key
                && record_end >= segment.start_1based
                && record_start <= segment.end_1based
        })
}

fn oligo_overlap(
    oligo: &PrimerVariantScreenOligo,
    record: &ParsedVcfRecord,
    critical_3prime_bases: usize,
    aliases: &BTreeMap<String, String>,
) -> Option<OligoOverlap> {
    let record_start = record.position_1based;
    let record_end = record_start.checked_add(record.reference_allele.len() - 1)?;
    let record_key = contig_key(&record.reference_name, aliases);
    let mut positions = Vec::new();
    let mut full_reference_match = None;
    for segment in &oligo.binding_segments {
        if contig_key(&segment.reference_name, aliases) != record_key
            || record_end < segment.start_1based
            || record_start > segment.end_1based
        {
            continue;
        }
        let overlap_start = record_start.max(segment.start_1based);
        let overlap_end = record_end.min(segment.end_1based);
        let mapped_start = map_genomic_to_oligo(segment, overlap_start);
        let mapped_end = map_genomic_to_oligo(segment, overlap_end);
        positions.push(mapped_start.min(mapped_end));
        positions.push(mapped_start.max(mapped_end));
        if record_start >= segment.start_1based && record_end <= segment.end_1based {
            let offset = record_start - segment.start_1based;
            let expected =
                &segment.reference_sequence_5_to_3[offset..offset + record.reference_allele.len()];
            full_reference_match = Some(expected == record.reference_allele);
        }
    }
    let first = positions.iter().min().copied()?;
    let last = positions.iter().max().copied()?;
    let distance = oligo.sequence_5_to_3.len().saturating_sub(last);
    Some(OligoOverlap {
        oligo_position_1based: first,
        oligo_end_position_1based: last,
        distance_from_3prime_end: distance,
        critical_three_prime: distance < critical_3prime_bases,
        reference_match: full_reference_match,
    })
}

fn amplicon_overlap(
    pair: &PairWork,
    record: &ParsedVcfRecord,
    aliases: &BTreeMap<String, String>,
) -> Option<Option<bool>> {
    let record_start = record.position_1based;
    let record_end = record_start.checked_add(record.reference_allele.len() - 1)?;
    let key = contig_key(&record.reference_name, aliases);
    let mut found = false;
    for segment in &pair.amplicon_segments {
        if contig_key(&segment.reference_name, aliases) != key
            || record_end < segment.start_1based
            || record_start > segment.end_1based
        {
            continue;
        }
        found = true;
        if record_start >= segment.start_1based && record_end <= segment.end_1based {
            let offset = record_start - segment.start_1based;
            let expected =
                &segment.reference_sequence_5_to_3[offset..offset + record.reference_allele.len()];
            return Some(Some(expected == record.reference_allele));
        }
    }
    found.then_some(None)
}

fn map_genomic_to_oligo(segment: &PrimerVariantBindingSegment, position_1based: usize) -> usize {
    match segment.strand {
        PrimerVariantBindingStrand::Plus => {
            segment.oligo_start_0based + (position_1based - segment.start_1based) + 1
        }
        PrimerVariantBindingStrand::Minus => {
            segment.oligo_start_0based + (segment.end_1based - position_1based) + 1
        }
    }
}

fn visit_vcf_records(
    path: &Path,
    mut visit_header: impl FnMut(&str),
    mut visit_record: impl FnMut(usize, ParsedVcfRecord) -> Result<(), String>,
) -> Result<(), String> {
    let mut reader = open_maybe_gz_reader(path)?;
    let mut line = String::new();
    let mut line_number = 0usize;
    let mut record_number = 0usize;
    while reader
        .read_line(&mut line)
        .map_err(|error| format!("Could not read VCF '{}': {error}", path.display()))?
        != 0
    {
        line_number += 1;
        let row = line.trim_end_matches(['\r', '\n']);
        if row.starts_with('#') {
            visit_header(row);
            line.clear();
            continue;
        }
        if row.trim().is_empty() {
            line.clear();
            continue;
        }
        let columns = row.split('\t').collect::<Vec<_>>();
        if columns.len() < 8 {
            return Err(format!(
                "VCF '{}' line {} has {} columns; expected at least 8",
                path.display(),
                line_number,
                columns.len()
            ));
        }
        let position_1based = columns[1].parse::<usize>().map_err(|error| {
            format!(
                "Invalid VCF position '{}' at '{}':{}: {error}",
                columns[1],
                path.display(),
                line_number
            )
        })?;
        if position_1based == 0 {
            return Err(format!(
                "VCF position must be >= 1 at '{}':{}",
                path.display(),
                line_number
            ));
        }
        let reference_allele = normalize_vcf_allele(columns[3], "REF", line_number)?;
        let alternate_alleles = columns[4]
            .split(',')
            .filter(|value| *value != "." && !value.is_empty())
            .map(|value| normalize_vcf_allele(value, "ALT", line_number))
            .collect::<Result<Vec<_>, _>>()?;
        if alternate_alleles.is_empty() {
            line.clear();
            continue;
        }
        record_number += 1;
        visit_record(
            record_number,
            ParsedVcfRecord {
                reference_name: columns[0].trim().to_string(),
                position_1based,
                id: columns[2].trim().to_string(),
                reference_allele,
                alternate_alleles,
                filter: columns[6].trim().to_string(),
                info: columns[7].trim().to_string(),
            },
        )?;
        line.clear();
    }
    Ok(())
}

fn open_maybe_gz_reader(path: &Path) -> Result<Box<dyn BufRead>, String> {
    let file = File::open(path)
        .map_err(|error| format!("Could not open VCF '{}': {error}", path.display()))?;
    if path
        .file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name.to_ascii_lowercase().ends_with(".gz"))
    {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn allele_frequencies(
    info: &str,
    field: Option<&str>,
    alternate_count: usize,
    warnings: &mut Vec<String>,
    record: &ParsedVcfRecord,
) -> Vec<Option<f64>> {
    let Some(field) = field else {
        return vec![None; alternate_count];
    };
    let raw = info.split(';').find_map(|token| {
        let (key, value) = token.split_once('=')?;
        key.eq_ignore_ascii_case(field).then_some(value)
    });
    let Some(raw) = raw else {
        return vec![None; alternate_count];
    };
    let parsed = raw
        .split(',')
        .map(|value| {
            if value == "." || value.trim().is_empty() {
                return None;
            }
            match value.parse::<f64>() {
                Ok(frequency) if (0.0..=1.0).contains(&frequency) => Some(frequency),
                _ => {
                    warnings.push(format!(
                        "Variant {}:{} has invalid {} value '{}'; frequency is unknown.",
                        record.reference_name, record.position_1based, field, value
                    ));
                    None
                }
            }
        })
        .collect::<Vec<_>>();
    (0..alternate_count)
        .map(|index| parsed.get(index).copied().flatten())
        .collect()
}

fn overlap_note(
    overlap_kind: PrimerVariantOverlapKind,
    variant_kind: PrimerVariantKind,
    frequency: Option<f64>,
    maximum_allowed_frequency: Option<f64>,
    critical_three_prime: bool,
    reference_match: Option<bool>,
) -> String {
    let frequency_text = match (frequency, maximum_allowed_frequency) {
        (Some(value), Some(threshold)) => {
            format!("allele frequency {value:.6} (allowed maximum {threshold:.6})")
        }
        (Some(value), None) => format!("allele frequency {value:.6}"),
        (None, _) => "allele frequency unknown (not treated as zero)".to_string(),
    };
    let reference_text = match reference_match {
        Some(true) => "reference allele verified",
        Some(false) => "reference allele mismatch",
        None => "reference allele could not be verified across the declared segment",
    };
    format!(
        "{:?} {:?} overlap; {}; {}; 3-prime-critical={critical_three_prime}.",
        overlap_kind, variant_kind, frequency_text, reference_text
    )
}

fn classify_variant(reference: &str, alternate: &str) -> PrimerVariantKind {
    if alternate.starts_with('<') || alternate == "*" || alternate.contains('[') {
        PrimerVariantKind::Complex
    } else {
        match (reference.len(), alternate.len()) {
            (1, 1) => PrimerVariantKind::Snv,
            (left, right) if left == right => PrimerVariantKind::Mnv,
            (left, right) if left < right => PrimerVariantKind::Insertion,
            (left, right) if left > right => PrimerVariantKind::Deletion,
            _ => PrimerVariantKind::Complex,
        }
    }
}

fn variant_id(record: &ParsedVcfRecord, alternate: &str, alternate_index: usize) -> String {
    if record.id.is_empty() || record.id == "." {
        format!(
            "{}:{}:{}>{}",
            record.reference_name, record.position_1based, record.reference_allele, alternate
        )
    } else if record.alternate_alleles.len() == 1 {
        record.id.clone()
    } else {
        format!("{}:alt{}", record.id, alternate_index + 1)
    }
}

fn parse_vcf_contig_header(header: &str) -> Option<String> {
    let body = header.strip_prefix("##contig=<")?.strip_suffix('>')?;
    body.split(',').find_map(|token| {
        let (key, value) = token.split_once('=')?;
        key.trim()
            .eq_ignore_ascii_case("ID")
            .then(|| value.trim().to_string())
    })
}

fn fingerprint_request(
    request: &PrimerVariantScreenRequest,
    source: &ResolvedVariantSource,
    content_sha256: &str,
) -> Result<String, String> {
    #[derive(Serialize)]
    struct Fingerprint<'a> {
        schema: &'a str,
        reference_assembly: &'a str,
        source_reference_assembly: &'a str,
        source_name: &'a str,
        source_release: &'a str,
        population: &'a str,
        retrieval_time: &'a str,
        allele_frequency_info_field: &'a Option<String>,
        contig_aliases: &'a BTreeMap<String, String>,
        manifest_sha256: &'a Option<String>,
        content_sha256: &'a str,
        candidates: &'a [PrimerVariantScreenCandidate],
        maximum_allowed_frequency: Option<f64>,
        critical_3prime_bases: usize,
        probe_overlap_policy: PrimerVariantProbeOverlapPolicy,
    }
    serde_json::to_vec(&Fingerprint {
        schema: PRIMER_VARIANT_SCREEN_REQUEST_SCHEMA,
        reference_assembly: request.reference_assembly.trim(),
        source_reference_assembly: &source.reference_assembly,
        source_name: &source.source_name,
        source_release: &source.source_release,
        population: &source.population,
        retrieval_time: &source.retrieval_time,
        allele_frequency_info_field: &source.allele_frequency_info_field,
        contig_aliases: &source.contig_aliases,
        manifest_sha256: &source.manifest_sha256,
        content_sha256,
        candidates: &request.candidates,
        maximum_allowed_frequency: request.maximum_allowed_frequency,
        critical_3prime_bases: request.critical_3prime_bases,
        probe_overlap_policy: request.probe_overlap_policy,
    })
    .map(|bytes| sha256_prefixed_bytes(&bytes))
    .map_err(|error| format!("Could not fingerprint primer variant request: {error}"))
}

fn fingerprint_geometry(
    forward: &PrimerVariantScreenOligo,
    reverse: &PrimerVariantScreenOligo,
    probe: Option<&PrimerVariantScreenOligo>,
    amplicon_segments: &[PrimerVariantAmpliconSegment],
) -> Result<String, String> {
    serde_json::to_vec(&(forward, reverse, probe, amplicon_segments))
        .map(|bytes| sha256_prefixed_bytes(&bytes))
        .map_err(|error| format!("Could not fingerprint primer binding geometry: {error}"))
}

fn contig_key(name: &str, aliases: &BTreeMap<String, String>) -> String {
    let mut current = name.trim().to_string();
    if let Some((_, canonical)) = aliases
        .iter()
        .find(|(alias, _)| alias.eq_ignore_ascii_case(&current))
    {
        current = canonical.clone();
    }
    let upper = current.to_ascii_uppercase();
    let stripped = upper.strip_prefix("CHR").unwrap_or(&upper);
    match stripped {
        "M" => "MT".to_string(),
        other => other.to_string(),
    }
}

fn assembly_matches(left: &str, right: &str) -> bool {
    left.trim().eq_ignore_ascii_case(right.trim())
}

fn normalize_dna(raw: &str, label: &str) -> Result<String, String> {
    let sequence = canonical_oligo_sequence(raw);
    if sequence.is_empty() {
        return Err(format!("{label} must not be empty"));
    }
    if let Some(base) = sequence
        .chars()
        .find(|base| !matches!(base, 'A' | 'C' | 'G' | 'T'))
    {
        return Err(format!("{label} contains unsupported base '{base}'"));
    }
    Ok(sequence)
}

fn normalize_vcf_allele(raw: &str, field: &str, line_number: usize) -> Result<String, String> {
    let allele = raw.trim().to_ascii_uppercase();
    if allele.is_empty() || allele == "." {
        return Err(format!("VCF line {line_number} has empty {field} allele"));
    }
    Ok(allele)
}

fn reverse_complement(sequence: &str) -> Result<String, String> {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => Ok('T'),
            b'C' => Ok('G'),
            b'G' => Ok('C'),
            b'T' => Ok('A'),
            _ => Err(format!("Cannot reverse-complement base '{}'", base as char)),
        })
        .collect()
}

fn validate_closed_range(start: usize, end: usize, label: &str) -> Result<(), String> {
    if start == 0 || end < start {
        return Err(format!(
            "{label} requires a valid 1-based closed range, got {start}-{end}"
        ));
    }
    Ok(())
}

fn required_text(raw: &str, field: &str) -> Result<String, String> {
    let value = raw.trim();
    if value.is_empty() {
        Err(format!("Primer variant screen requires non-empty {field}"))
    } else {
        Ok(value.to_string())
    }
}

fn normalize_optional_text(raw: Option<&str>) -> Option<String> {
    raw.map(str::trim)
        .filter(|value| !value.is_empty())
        .map(str::to_string)
}

fn require_matching_override(raw: &str, manifest: &str, field: &str) -> Result<(), String> {
    let raw = raw.trim();
    if !raw.is_empty() && !raw.eq_ignore_ascii_case(manifest.trim()) {
        Err(format!(
            "Request {field} '{}' conflicts with manifest value '{}'",
            raw, manifest
        ))
    } else {
        Ok(())
    }
}

fn require_matching_optional_override(
    raw: Option<&str>,
    manifest: Option<&str>,
    field: &str,
) -> Result<(), String> {
    let raw = normalize_optional_text(raw);
    let manifest = normalize_optional_text(manifest);
    if raw.is_some() && raw != manifest {
        Err(format!(
            "Request {field} {:?} conflicts with manifest value {:?}",
            raw, manifest
        ))
    } else {
        Ok(())
    }
}

fn normalize_sha256(raw: &str) -> String {
    let value = raw.trim().to_ascii_lowercase();
    if value.starts_with("sha256:") {
        value
    } else {
        format!("sha256:{value}")
    }
}

fn ensure_vcf_exists(path: &Path) -> Result<(), String> {
    if path.is_file() {
        Ok(())
    } else {
        Err(format!("Variant VCF '{}' does not exist", path.display()))
    }
}

fn short_digest_id(prefix: &str, digest: &str) -> String {
    let hex = digest.strip_prefix("sha256:").unwrap_or(digest);
    format!("{prefix}_{}", hex.chars().take(16).collect::<String>())
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{Compression, write::GzEncoder};
    use std::io::Write;
    use tempfile::{TempDir, tempdir};

    const FORWARD: &str = "AACCGGTTAACC";
    const REVERSE: &str = "GGAACCTTGGAA";
    const REVERSE_REFERENCE: &str = "TTCCAAGGTTCC";
    const PROBE: &str = "ACGTACGT";

    fn write_vcf(dir: &TempDir, rows: &str) -> PathBuf {
        let path = dir.path().join("variants.vcf");
        std::fs::write(
            &path,
            format!(
                "##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{rows}"
            ),
        )
        .expect("write synthetic VCF");
        path
    }

    fn gzip_file(path: &Path) -> PathBuf {
        let gzip_path = path.with_extension("vcf.gz");
        let input = std::fs::read(path).expect("read VCF for gzip");
        let file = File::create(&gzip_path).expect("create gzip VCF");
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(&input).expect("write gzip VCF");
        encoder.finish().expect("finish gzip VCF");
        gzip_path
    }

    fn amplicon_reference() -> String {
        let mut sequence = vec![b'A'; 112];
        sequence[0..12].copy_from_slice(FORWARD.as_bytes());
        sequence[50..58].copy_from_slice(PROBE.as_bytes());
        sequence[100..112].copy_from_slice(REVERSE_REFERENCE.as_bytes());
        String::from_utf8(sequence).expect("ASCII amplicon")
    }

    fn standard_candidate(candidate_id: &str) -> PrimerVariantScreenCandidate {
        PrimerVariantScreenCandidate {
            candidate_id: candidate_id.to_string(),
            source: PrimerVariantCandidateSource {
                candidate_id: candidate_id.to_string(),
                source_kind: crate::engine::PrimerVariantCandidateSourceKind::DeNovo,
                source_id: candidate_id.to_string(),
                provider: "GENtle".to_string(),
                content_sha256: format!("sha256:{candidate_id}"),
                ..Default::default()
            },
            forward: PrimerVariantScreenOligo {
                sequence_5_to_3: FORWARD.to_string(),
                binding_segments: vec![PrimerVariantBindingSegment {
                    reference_name: "chr1".to_string(),
                    start_1based: 100,
                    end_1based: 111,
                    strand: PrimerVariantBindingStrand::Plus,
                    oligo_start_0based: 0,
                    oligo_end_0based_exclusive: 12,
                    reference_sequence_5_to_3: FORWARD.to_string(),
                }],
            },
            reverse: PrimerVariantScreenOligo {
                sequence_5_to_3: REVERSE.to_string(),
                binding_segments: vec![PrimerVariantBindingSegment {
                    reference_name: "chr1".to_string(),
                    start_1based: 200,
                    end_1based: 211,
                    strand: PrimerVariantBindingStrand::Minus,
                    oligo_start_0based: 0,
                    oligo_end_0based_exclusive: 12,
                    reference_sequence_5_to_3: REVERSE_REFERENCE.to_string(),
                }],
            },
            probe: Some(PrimerVariantScreenOligo {
                sequence_5_to_3: PROBE.to_string(),
                binding_segments: vec![PrimerVariantBindingSegment {
                    reference_name: "chr1".to_string(),
                    start_1based: 150,
                    end_1based: 157,
                    strand: PrimerVariantBindingStrand::Plus,
                    oligo_start_0based: 0,
                    oligo_end_0based_exclusive: 8,
                    reference_sequence_5_to_3: PROBE.to_string(),
                }],
            }),
            amplicon_segments: vec![PrimerVariantAmpliconSegment {
                reference_name: "chr1".to_string(),
                start_1based: 100,
                end_1based: 211,
                reference_sequence_5_to_3: amplicon_reference(),
            }],
            ..Default::default()
        }
    }

    fn request(
        path: &Path,
        candidates: Vec<PrimerVariantScreenCandidate>,
    ) -> PrimerVariantScreenRequest {
        PrimerVariantScreenRequest {
            reference_assembly: "GRCh38".to_string(),
            source: PrimerVariantSourceInput {
                vcf_path: Some(path.display().to_string()),
                reference_assembly: "GRCh38".to_string(),
                source_name: "synthetic variants".to_string(),
                source_release: "1".to_string(),
                population: "synthetic".to_string(),
                retrieval_time: "deterministic-fixture".to_string(),
                allele_frequency_info_field: Some("AF".to_string()),
                ..Default::default()
            },
            candidates,
            maximum_allowed_frequency: Some(0.01),
            ..Default::default()
        }
    }

    #[test]
    fn maps_plus_and_minus_terminal_variants_strand_aware() {
        let dir = tempdir().expect("temp dir");
        let vcf = write_vcf(
            &dir,
            "chr1\t111\trsPlus3p\tC\tT\t.\tPASS\tAF=0.2\nchr1\t200\trsMinus3p\tT\tC\t.\tPASS\tAF=0.3\n",
        );
        let report = screen_primer_variants(request(&vcf, vec![standard_candidate("denovo")]))
            .expect("variant screen");
        let evidence = &report.evidence_reports[0];
        assert_eq!(
            evidence.status,
            PrimerVariantEvidenceStatus::VariantDetected
        );
        assert_eq!(evidence.overlaps.len(), 2);
        for overlap in &evidence.overlaps {
            assert_eq!(overlap.oligo_position_1based, Some(12));
            assert_eq!(overlap.distance_from_3prime_end, Some(0));
            assert!(overlap.critical_three_prime);
            assert_eq!(overlap.reference_match, Some(true));
            assert!(overlap.relevant_under_policy);
        }
        assert_eq!(evidence.overlaps[0].role, "forward");
        assert_eq!(evidence.overlaps[1].role, "reverse");
    }

    #[test]
    fn separates_rare_primer_probe_and_amplicon_only_variation() {
        let dir = tempdir().expect("temp dir");
        let vcf = write_vcf(
            &dir,
            concat!(
                "chr1\t104\trsRareInternal\tG\tA\t.\tPASS\tAF=0.001\n",
                "chr1\t153\trsProbeUnknown\tT\tC\t.\tPASS\t.\n",
                "chr1\t170\trsAmplicon\tA\tG\t.\tPASS\tAF=0.4\n"
            ),
        );
        let report = screen_primer_variants(request(&vcf, vec![standard_candidate("mixed")]))
            .expect("variant screen");
        let rows = &report.evidence_reports[0].overlaps;
        let rare = rows
            .iter()
            .find(|row| row.variant_id == "rsRareInternal")
            .expect("rare primer row");
        assert_eq!(rare.overlap_kind, PrimerVariantOverlapKind::Primer);
        assert_eq!(rare.allele_frequency, Some(0.001));
        assert!(!rare.relevant_under_policy);
        let probe = rows
            .iter()
            .find(|row| row.variant_id == "rsProbeUnknown")
            .expect("probe row");
        assert_eq!(probe.overlap_kind, PrimerVariantOverlapKind::Probe);
        assert_eq!(probe.allele_frequency, None);
        assert!(
            probe.relevant_under_policy,
            "unknown AF must not become zero"
        );
        let amplicon = rows
            .iter()
            .find(|row| row.variant_id == "rsAmplicon")
            .expect("amplicon row");
        assert_eq!(
            amplicon.overlap_kind,
            PrimerVariantOverlapKind::AmpliconOnly
        );
        assert!(!amplicon.relevant_under_policy);
    }

    #[test]
    fn probe_policy_report_only_and_ignore_remain_distinct_from_amplicon_only() {
        let dir = tempdir().expect("temp dir");
        let vcf = write_vcf(&dir, "chr1\t153\trsProbe\tT\tC\t.\tPASS\tAF=0.2\n");
        let mut report_only_request = request(&vcf, vec![standard_candidate("report-only")]);
        report_only_request.probe_overlap_policy = PrimerVariantProbeOverlapPolicy::ReportOnly;
        let report_only =
            screen_primer_variants(report_only_request).expect("report-only probe screen");
        assert_eq!(report_only.evidence_reports[0].overlaps.len(), 1);
        assert_eq!(
            report_only.evidence_reports[0].overlaps[0].overlap_kind,
            PrimerVariantOverlapKind::Probe
        );
        assert!(!report_only.evidence_reports[0].overlaps[0].relevant_under_policy);
        assert_eq!(
            report_only.evidence_reports[0].status,
            PrimerVariantEvidenceStatus::EvaluatedNoRelevantVariant
        );

        let mut ignored_request = request(&vcf, vec![standard_candidate("ignored")]);
        ignored_request.probe_overlap_policy = PrimerVariantProbeOverlapPolicy::Ignore;
        let ignored = screen_primer_variants(ignored_request).expect("ignored probe screen");
        assert!(ignored.evidence_reports[0].overlaps.is_empty());
        assert_eq!(
            ignored.evidence_reports[0].status,
            PrimerVariantEvidenceStatus::EvaluatedNoRelevantVariant
        );
    }

    #[test]
    fn junction_binding_segments_map_to_one_oligo_and_retain_duplicate_sources() {
        let dir = tempdir().expect("temp dir");
        let vcf = gzip_file(&write_vcf(
            &dir,
            "chr1\t202\trsJunctionExon\tC\tT\t.\tPASS\tAF=0.2\n",
        ));
        let mut first = standard_candidate("primerbank-row");
        first.forward = PrimerVariantScreenOligo {
            sequence_5_to_3: "AAAACCCC".to_string(),
            binding_segments: vec![
                PrimerVariantBindingSegment {
                    reference_name: "chr1".to_string(),
                    start_1based: 100,
                    end_1based: 103,
                    strand: PrimerVariantBindingStrand::Plus,
                    oligo_start_0based: 0,
                    oligo_end_0based_exclusive: 4,
                    reference_sequence_5_to_3: "AAAA".to_string(),
                },
                PrimerVariantBindingSegment {
                    reference_name: "chr1".to_string(),
                    start_1based: 200,
                    end_1based: 203,
                    strand: PrimerVariantBindingStrand::Plus,
                    oligo_start_0based: 4,
                    oligo_end_0based_exclusive: 8,
                    reference_sequence_5_to_3: "CCCC".to_string(),
                },
            ],
        };
        first.amplicon_segments.clear();
        first.source.source_kind = crate::engine::PrimerVariantCandidateSourceKind::PrimerBank;
        let mut second = first.clone();
        second.candidate_id = "vendor-row".to_string();
        second.source.candidate_id = second.candidate_id.clone();
        second.source.source_kind =
            crate::engine::PrimerVariantCandidateSourceKind::CommercialCatalogue;
        second.source.source_id = "CAT-1".to_string();

        let manifest_path = dir.path().join("resource.json");
        let manifest = PrimerVariantResourceManifest {
            schema: PRIMER_VARIANT_RESOURCE_MANIFEST_SCHEMA.to_string(),
            vcf_path: vcf
                .file_name()
                .expect("VCF file name")
                .to_string_lossy()
                .to_string(),
            reference_assembly: "GRCh38".to_string(),
            source_name: "synthetic variants".to_string(),
            source_release: "1".to_string(),
            population: "synthetic".to_string(),
            retrieval_time: "deterministic-fixture".to_string(),
            allele_frequency_info_field: Some("AF".to_string()),
            ..Default::default()
        };
        std::fs::write(
            &manifest_path,
            serde_json::to_vec_pretty(&manifest).expect("serialize manifest"),
        )
        .expect("write manifest");
        let mut request = request(&vcf, vec![first, second]);
        request.source = PrimerVariantSourceInput {
            manifest_path: Some(manifest_path.display().to_string()),
            ..Default::default()
        };
        let report = screen_primer_variants(request).expect("manifest gzip screen");
        assert_eq!(report.candidate_count, 2);
        assert_eq!(report.unique_pair_count, 1);
        let evidence = &report.evidence_reports[0];
        assert_eq!(evidence.candidate_sources.len(), 2);
        assert_eq!(evidence.overlaps[0].oligo_position_1based, Some(7));
        assert_eq!(evidence.overlaps[0].reference_match, Some(true));
    }

    #[test]
    fn assembly_and_reference_mismatches_fail_closed_and_indels_warn() {
        let dir = tempdir().expect("temp dir");
        let mismatch_vcf = write_vcf(&dir, "chr1\t111\trsWrongRef\tA\tT\t.\tPASS\tAF=0.2\n");
        let mismatch = screen_primer_variants(request(
            &mismatch_vcf,
            vec![standard_candidate("wrong-ref")],
        ))
        .expect("reference mismatch is reportable");
        assert_eq!(
            mismatch.evidence_reports[0].status,
            PrimerVariantEvidenceStatus::IncompatibleReference
        );
        assert_eq!(
            mismatch.evidence_reports[0].overlaps[0].reference_match,
            Some(false)
        );

        let mut wrong_assembly_request =
            request(&mismatch_vcf, vec![standard_candidate("wrong-assembly")]);
        wrong_assembly_request.source.reference_assembly = "GRCh37".to_string();
        let wrong_assembly = screen_primer_variants(wrong_assembly_request)
            .expect("assembly mismatch is reportable");
        assert_eq!(wrong_assembly.vcf_record_count, 0);
        assert_eq!(
            wrong_assembly.evidence_reports[0].status,
            PrimerVariantEvidenceStatus::IncompatibleReference
        );

        let wrong_contig_vcf = dir.path().join("wrong-contig.vcf");
        std::fs::write(
            &wrong_contig_vcf,
            "##fileformat=VCFv4.2\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr2\t111\trsOtherContig\tC\tT\t.\tPASS\tAF=0.2\n",
        )
        .expect("write wrong-contig VCF");
        let wrong_contig = screen_primer_variants(request(
            &wrong_contig_vcf,
            vec![standard_candidate("wrong-contig")],
        ))
        .expect("contig mismatch is reportable");
        assert_eq!(
            wrong_contig.evidence_reports[0].status,
            PrimerVariantEvidenceStatus::IncompatibleReference
        );
        assert!(
            wrong_contig.evidence_reports[0]
                .warnings
                .iter()
                .any(|warning| warning.contains("observed records"))
        );

        let indel_vcf = write_vcf(&dir, "chr1\t105\trsInsertion\tG\tGA\t.\tPASS\tAF=0.2\n");
        let indel = screen_primer_variants(request(&indel_vcf, vec![standard_candidate("indel")]))
            .expect("indel screen");
        assert_eq!(
            indel.evidence_reports[0].overlaps[0].variant_kind,
            PrimerVariantKind::Insertion
        );
        assert!(
            indel.evidence_reports[0]
                .warnings
                .iter()
                .any(|warning| warning.contains("conservatively"))
        );
    }

    #[test]
    fn old_evidence_payload_defaults_new_additive_fields() {
        let old = serde_json::json!({
            "schema": "gentle.primer_variant_evidence.v1",
            "evidence_id": "old",
            "pair_id": "pair_sha256_old",
            "reference_assembly": "GRCh38",
            "source_name": "old source",
            "source_release": "1",
            "population": "all",
            "normalization_method": "old",
            "retrieval_time": "old",
            "content_sha256": "sha256:old",
            "status": "variant_detected",
            "overlaps": [{
                "oligo_id": "oligo_sha256_old",
                "role": "forward",
                "variant_id": "rsOld",
                "reference_name": "chr1",
                "position_1based": 10,
                "reference_allele": "A",
                "alternate_allele": "G",
                "relevant_under_policy": true,
                "note": "legacy"
            }]
        });
        let report: PrimerVariantEvidenceReport =
            serde_json::from_value(old).expect("old evidence remains readable");
        assert_eq!(report.critical_3prime_bases, 5);
        assert!(report.candidate_sources.is_empty());
        assert_eq!(
            report.overlaps[0].overlap_kind,
            PrimerVariantOverlapKind::Primer
        );
        assert_eq!(report.overlaps[0].reference_end_1based, 0);
    }
}
