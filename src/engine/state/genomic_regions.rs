//! Engine-owned persistence, projection, and interchange for genomic ROIs.

use super::*;
use gentle_protocol as gp;

const BED_COLUMNS: [&str; 6] = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"];

fn region_error(code: ErrorCode, message: impl Into<String>) -> EngineError {
    EngineError {
        code,
        message: message.into(),
        cause_chain: vec![],
    }
}

fn nonempty(value: &str, field: &str) -> Result<(), EngineError> {
    if value.trim().is_empty() || value.chars().any(char::is_control) {
        return Err(region_error(
            ErrorCode::InvalidInput,
            format!("{field} must be non-empty and contain no control characters"),
        ));
    }
    Ok(())
}

fn valid_identifier(value: &str, field: &str) -> Result<(), EngineError> {
    nonempty(value, field)?;
    if value.len() > 200
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.' | b':'))
    {
        return Err(region_error(
            ErrorCode::InvalidInput,
            format!("{field} must use at most 200 ASCII identifier characters"),
        ));
    }
    Ok(())
}

fn validate_reference(reference: &gp::GenomicRegionReference) -> Result<(), EngineError> {
    nonempty(&reference.assembly_name, "reference.assembly_name")?;
    nonempty(&reference.contig_name, "reference.contig_name")?;
    if reference
        .assembly_accession
        .as_deref()
        .is_some_and(|value| value.trim().is_empty())
        || reference
            .contig_accession
            .as_deref()
            .is_some_and(|value| value.trim().is_empty())
    {
        return Err(region_error(
            ErrorCode::InvalidInput,
            "reference accessions must be absent or non-empty",
        ));
    }
    Ok(())
}

fn validate_interval(interval: &gp::GenomicRegionInterval) -> Result<(), EngineError> {
    validate_reference(&interval.reference)?;
    if interval.start_0based >= interval.end_0based_exclusive {
        return Err(region_error(
            ErrorCode::InvalidInput,
            format!(
                "invalid 0-based half-open interval {}..{}; require start < end",
                interval.start_0based, interval.end_0based_exclusive
            ),
        ));
    }
    Ok(())
}

fn normalize_display_color(value: Option<String>) -> Result<Option<String>, EngineError> {
    let Some(value) = value else {
        return Ok(None);
    };
    let value = value.trim();
    if value.len() != 7
        || !value.starts_with('#')
        || !value[1..].bytes().all(|byte| byte.is_ascii_hexdigit())
    {
        return Err(region_error(
            ErrorCode::InvalidInput,
            "display_color_hex must be '#RRGGBB'",
        ));
    }
    Ok(Some(value.to_ascii_uppercase()))
}

fn validate_sha256(value: &str, field: &str) -> Result<(), EngineError> {
    let Some(hex) = value.strip_prefix("sha256:") else {
        return Err(region_error(
            ErrorCode::InvalidInput,
            format!("{field} must use the 'sha256:' prefix"),
        ));
    };
    if hex.len() != 64 || !hex.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        return Err(region_error(
            ErrorCode::InvalidInput,
            format!("{field} must contain exactly 64 hexadecimal SHA-256 characters"),
        ));
    }
    Ok(())
}

fn validate_local_projection(
    interval: &gp::GenomicRegionInterval,
    projection: &gp::GenomicRegionLocalProjection,
) -> Result<(), EngineError> {
    nonempty(&projection.seq_id, "local_projection.seq_id")?;
    validate_sha256(
        &projection.sequence_sha256,
        "local_projection.sequence_sha256",
    )?;
    nonempty(
        &projection.source_genome_id,
        "local_projection.source_genome_id",
    )?;
    if projection.anchor_start_1based == 0
        || projection.anchor_start_1based > projection.anchor_end_1based
    {
        return Err(region_error(
            ErrorCode::InvalidInput,
            "local projection requires a non-empty 1-based inclusive genome anchor",
        ));
    }
    if projection.anchor_strand == gp::GenomicRegionStrand::Unstranded {
        return Err(region_error(
            ErrorCode::InvalidInput,
            "local projection genome anchor must have an explicit plus or minus strand",
        ));
    }
    let anchor_len = projection
        .anchor_end_1based
        .checked_sub(projection.anchor_start_1based)
        .and_then(|span| span.checked_add(1))
        .ok_or_else(|| region_error(ErrorCode::InvalidInput, "invalid local projection anchor"))?;
    if projection.local_start_0based >= projection.local_end_0based_exclusive
        || projection.local_end_0based_exclusive > anchor_len
    {
        return Err(region_error(
            ErrorCode::InvalidInput,
            "local projection interval falls outside its exact genome anchor",
        ));
    }
    let assembly_matches = projection.source_genome_id == interval.reference.assembly_name
        || interval.reference.assembly_accession.as_deref()
            == Some(projection.source_genome_id.as_str());
    if !assembly_matches {
        return Err(region_error(
            ErrorCode::InvalidInput,
            "local projection genome id does not match the canonical interval assembly",
        ));
    }
    let anchor_start_0based = projection.anchor_start_1based - 1;
    let anchor_end_0based_exclusive = projection.anchor_end_1based;
    let (expected_start, expected_end, expected_strand) =
        if projection.anchor_strand == gp::GenomicRegionStrand::Minus {
            (
                anchor_end_0based_exclusive
                    .checked_sub(projection.local_end_0based_exclusive)
                    .ok_or_else(|| {
                        region_error(
                            ErrorCode::InvalidInput,
                            "local projection underflows anchor",
                        )
                    })?,
                anchor_end_0based_exclusive
                    .checked_sub(projection.local_start_0based)
                    .ok_or_else(|| {
                        region_error(
                            ErrorCode::InvalidInput,
                            "local projection underflows anchor",
                        )
                    })?,
                flip_strand(projection.local_strand),
            )
        } else {
            (
                anchor_start_0based
                    .checked_add(projection.local_start_0based)
                    .ok_or_else(|| {
                        region_error(ErrorCode::InvalidInput, "local projection overflows anchor")
                    })?,
                anchor_start_0based
                    .checked_add(projection.local_end_0based_exclusive)
                    .ok_or_else(|| {
                        region_error(ErrorCode::InvalidInput, "local projection overflows anchor")
                    })?,
                projection.local_strand,
            )
        };
    if (expected_start, expected_end, expected_strand)
        != (
            interval.start_0based,
            interval.end_0based_exclusive,
            interval.strand,
        )
    {
        return Err(region_error(
            ErrorCode::InvalidInput,
            "local projection geometry or strand does not reproduce the canonical genomic interval",
        ));
    }
    Ok(())
}

fn serialize_for_digest<T: Serialize>(value: &T, label: &str) -> Result<String, EngineError> {
    serde_json::to_string(value).map_err(|error| {
        region_error(
            ErrorCode::Internal,
            format!("could not serialize {label} for deterministic digest: {error}"),
        )
    })
}

fn sort_region_evidence(evidence: &mut [gp::GenomicRegionEvidenceReference]) {
    evidence.sort_by(|left, right| {
        (&left.evidence_id, &left.source_kind, &left.source_id).cmp(&(
            &right.evidence_id,
            &right.source_kind,
            &right.source_id,
        ))
    });
}

fn recompute_region_digests(region: &mut gp::GenomicRegionOfInterest) -> Result<(), EngineError> {
    sort_region_evidence(&mut region.evidence);
    if let Some(derivation) = region.derivation.as_mut() {
        derivation.parents.sort_by(|left, right| {
            (&left.region_id, &left.identity_sha256)
                .cmp(&(&right.region_id, &right.identity_sha256))
        });
    }
    let identity_payload = serde_json::json!({
        "schema": region.schema,
        "interval": region.interval,
        "local_projection": region.local_projection.as_ref().map(|projection| serde_json::json!({
            "seq_id": projection.seq_id,
            "sequence_sha256": projection.sequence_sha256,
            "source_genome_id": projection.source_genome_id,
            "anchor_start_1based": projection.anchor_start_1based,
            "anchor_end_1based": projection.anchor_end_1based,
            "anchor_strand": projection.anchor_strand,
            "local_start_0based": projection.local_start_0based,
            "local_end_0based_exclusive": projection.local_end_0based_exclusive,
            "local_strand": projection.local_strand,
        })),
        "purpose": region.purpose,
        "selection_method": region.selection_method,
        "evidence": region.evidence,
        "derivation": region.derivation,
    });
    region.identity_sha256 = sha256_prefixed_str(&serialize_for_digest(
        &identity_payload,
        "genomic-region identity",
    )?);
    if region.region_id.trim().is_empty() {
        region.region_id = format!(
            "region_{}",
            region
                .identity_sha256
                .strip_prefix("sha256:")
                .unwrap_or(&region.identity_sha256)
                .chars()
                .take(16)
                .collect::<String>()
        );
    }
    let mut content = region.clone();
    content.content_sha256.clear();
    region.content_sha256 =
        sha256_prefixed_str(&serialize_for_digest(&content, "genomic-region content")?);
    Ok(())
}

fn recompute_set_digest(set: &mut gp::GenomicRegionSet) -> Result<(), EngineError> {
    set.regions
        .sort_by(|left, right| left.region_id.cmp(&right.region_id));
    let mut content = set.clone();
    content.content_sha256.clear();
    set.content_sha256 =
        sha256_prefixed_str(&serialize_for_digest(&content, "genomic-region set")?);
    Ok(())
}

fn validate_region_digest(region: &gp::GenomicRegionOfInterest) -> Result<(), EngineError> {
    validate_interval(&region.interval)?;
    valid_identifier(&region.region_id, "region_id")?;
    if region.schema != gp::GENOMIC_REGION_OF_INTEREST_SCHEMA {
        return Err(region_error(
            ErrorCode::InvalidInput,
            format!(
                "unsupported genomic ROI schema '{}'; expected '{}'",
                region.schema,
                gp::GENOMIC_REGION_OF_INTEREST_SCHEMA
            ),
        ));
    }
    if let Some(projection) = region.local_projection.as_ref() {
        validate_local_projection(&region.interval, projection)?;
    }
    let mut evidence_ids = BTreeSet::new();
    for evidence in &region.evidence {
        nonempty(&evidence.evidence_id, "evidence.evidence_id")?;
        nonempty(&evidence.source_kind, "evidence.source_kind")?;
        nonempty(&evidence.source_id, "evidence.source_id")?;
        if !evidence_ids.insert(evidence.evidence_id.as_str()) {
            return Err(region_error(
                ErrorCode::InvalidInput,
                format!("duplicate evidence_id '{}'", evidence.evidence_id),
            ));
        }
        if let Some(digest) = evidence.source_sha256.as_deref() {
            validate_sha256(digest, "evidence.source_sha256")?;
        }
    }
    if let Some(derivation) = region.derivation.as_ref() {
        if derivation.parents.is_empty() {
            return Err(region_error(
                ErrorCode::InvalidInput,
                "derived genomic region must retain at least one parent binding",
            ));
        }
        let mut parent_ids = BTreeSet::new();
        for parent in &derivation.parents {
            valid_identifier(&parent.region_id, "derivation.parent.region_id")?;
            validate_sha256(&parent.identity_sha256, "derivation.parent.identity_sha256")?;
            validate_sha256(&parent.content_sha256, "derivation.parent.content_sha256")?;
            if !parent_ids.insert(parent.region_id.as_str()) {
                return Err(region_error(
                    ErrorCode::InvalidInput,
                    format!("duplicate derivation parent '{}'", parent.region_id),
                ));
            }
        }
    }
    let mut expected = region.clone();
    recompute_region_digests(&mut expected)?;
    if expected.identity_sha256 != region.identity_sha256
        || expected.content_sha256 != region.content_sha256
    {
        return Err(region_error(
            ErrorCode::InvalidInput,
            format!("digest mismatch for genomic region '{}'", region.region_id),
        ));
    }
    Ok(())
}

fn validate_set_digest(set: &gp::GenomicRegionSet) -> Result<(), EngineError> {
    valid_identifier(&set.set_id, "set_id")?;
    if set.schema != gp::GENOMIC_REGION_SET_SCHEMA {
        return Err(region_error(
            ErrorCode::InvalidInput,
            format!(
                "unsupported genomic region-set schema '{}'; expected '{}'",
                set.schema,
                gp::GENOMIC_REGION_SET_SCHEMA
            ),
        ));
    }
    let mut ids = BTreeSet::new();
    let mut set_reference: Option<&gp::GenomicRegionReference> = None;
    for region in &set.regions {
        validate_region_digest(region)?;
        if let Some(reference) = set_reference {
            let candidate = &region.interval.reference;
            let incompatible_assembly = reference.assembly_name != candidate.assembly_name
                || matches!(
                    (&reference.assembly_accession, &candidate.assembly_accession),
                    (Some(left), Some(right)) if left != right
                );
            let incompatible_species = matches!(
                (&reference.species_scientific_name, &candidate.species_scientific_name),
                (Some(left), Some(right)) if left != right
            ) || matches!(
                (reference.taxon_id, candidate.taxon_id),
                (Some(left), Some(right)) if left != right
            );
            if incompatible_assembly || incompatible_species {
                return Err(region_error(
                    ErrorCode::InvalidInput,
                    format!(
                        "region set '{}' mixes incompatible species or assemblies; split the regions into separate sets instead",
                        set.set_id
                    ),
                ));
            }
        } else {
            set_reference = Some(&region.interval.reference);
        }
        if !ids.insert(region.region_id.as_str()) {
            return Err(region_error(
                ErrorCode::InvalidInput,
                format!("duplicate region_id '{}'", region.region_id),
            ));
        }
    }
    let mut expected = set.clone();
    recompute_set_digest(&mut expected)?;
    if expected.content_sha256 != set.content_sha256 {
        return Err(region_error(
            ErrorCode::InvalidInput,
            format!("digest mismatch for genomic region set '{}'", set.set_id),
        ));
    }
    Ok(())
}

fn references_match(left: &gp::GenomicRegionReference, right: &gp::GenomicRegionReference) -> bool {
    let assembly_matches = left.assembly_name == right.assembly_name
        && match (&left.assembly_accession, &right.assembly_accession) {
            (Some(left), Some(right)) => left == right,
            _ => true,
        };
    let left_contigs = std::iter::once(left.contig_name.as_str())
        .chain(left.contig_accession.as_deref())
        .chain(left.contig_aliases.iter().map(String::as_str))
        .collect::<BTreeSet<_>>();
    let right_contigs = std::iter::once(right.contig_name.as_str())
        .chain(right.contig_accession.as_deref())
        .chain(right.contig_aliases.iter().map(String::as_str))
        .collect::<BTreeSet<_>>();
    assembly_matches && !left_contigs.is_disjoint(&right_contigs)
}

fn strand_from_char(strand: Option<char>) -> gp::GenomicRegionStrand {
    match strand {
        Some('+') => gp::GenomicRegionStrand::Plus,
        Some('-') => gp::GenomicRegionStrand::Minus,
        _ => gp::GenomicRegionStrand::Unstranded,
    }
}

fn flip_strand(strand: gp::GenomicRegionStrand) -> gp::GenomicRegionStrand {
    match strand {
        gp::GenomicRegionStrand::Plus => gp::GenomicRegionStrand::Minus,
        gp::GenomicRegionStrand::Minus => gp::GenomicRegionStrand::Plus,
        gp::GenomicRegionStrand::Unstranded => gp::GenomicRegionStrand::Unstranded,
    }
}

fn human_copy(region: &gp::GenomicRegionOfInterest) -> String {
    let reference = &region.interval.reference;
    let accession = reference
        .assembly_accession
        .as_deref()
        .map(|value| format!(" ({value})"))
        .unwrap_or_default();
    let label = region.label.as_deref().unwrap_or(&region.region_id);
    let sources = region
        .evidence
        .iter()
        .map(|evidence| format!("{}:{}", evidence.source_kind, evidence.source_id))
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let source_summary = if sources.is_empty() {
        "none".to_string()
    } else {
        sources.join(",")
    };
    format!(
        "{} [{}] | {}{}, {}:{}-{}, 1-based inclusive | strand {} | method {} | sources {}",
        label,
        region.region_id,
        reference.assembly_name,
        accession,
        reference.contig_name,
        region.interval.start_0based + 1,
        region.interval.end_0based_exclusive,
        region.interval.strand.human_value(),
        region.selection_method.as_str(),
        source_summary
    )
}

fn bed_label(region: &gp::GenomicRegionOfInterest) -> String {
    let source = region.label.as_deref().unwrap_or(&region.region_id);
    let sanitized = source
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.' | ':') {
                ch
            } else {
                '_'
            }
        })
        .take(200)
        .collect::<String>();
    if sanitized.is_empty() {
        region.region_id.clone()
    } else {
        sanitized
    }
}

fn bed_row(region: &gp::GenomicRegionOfInterest) -> String {
    format!(
        "{}\t{}\t{}\t{}\t0\t{}",
        region.interval.reference.contig_name,
        region.interval.start_0based,
        region.interval.end_0based_exclusive,
        bed_label(region),
        region.interval.strand.bed_value()
    )
}

fn bed_payload(set: &gp::GenomicRegionSet) -> String {
    let mut payload = set
        .regions
        .iter()
        .map(bed_row)
        .collect::<Vec<_>>()
        .join("\n");
    if !payload.is_empty() {
        payload.push('\n');
    }
    payload
}

fn report_for_region(
    action: &str,
    set: Option<gp::GenomicRegionSet>,
    region: Option<gp::GenomicRegionOfInterest>,
) -> Result<gp::GenomicRegionOperationReport, EngineError> {
    let (human_copy, bed_row, canonical_roi_json) = if let Some(region) = region.as_ref() {
        (
            Some(human_copy(region)),
            Some(bed_row(region)),
            Some(serde_json::to_string_pretty(region).map_err(|error| {
                region_error(
                    ErrorCode::Internal,
                    format!("could not serialize canonical ROI JSON: {error}"),
                )
            })?),
        )
    } else {
        (None, None, None)
    };
    Ok(gp::GenomicRegionOperationReport {
        schema: gp::GENOMIC_REGION_OPERATION_REPORT_SCHEMA.to_string(),
        action: action.to_string(),
        set,
        region,
        sets: vec![],
        human_copy,
        bed_row,
        canonical_roi_json,
        written_artifacts: vec![],
        warnings: vec![],
    })
}

impl GentleEngine {
    fn genomic_region_store(&self) -> Result<gp::GenomicRegionStore, EngineError> {
        let Some(value) = self.state.metadata.get(GENOMIC_REGION_SETS_METADATA_KEY) else {
            return Ok(gp::GenomicRegionStore {
                schema: gp::GENOMIC_REGION_STORE_SCHEMA.to_string(),
                sets: vec![],
            });
        };
        let store: gp::GenomicRegionStore =
            serde_json::from_value(value.clone()).map_err(|error| {
                region_error(
                    ErrorCode::InvalidInput,
                    format!("could not decode persisted genomic-region store: {error}"),
                )
            })?;
        if store.schema != gp::GENOMIC_REGION_STORE_SCHEMA {
            return Err(region_error(
                ErrorCode::InvalidInput,
                format!(
                    "unsupported persisted genomic-region store schema '{}'",
                    store.schema
                ),
            ));
        }
        for set in &store.sets {
            validate_set_digest(set)?;
        }
        Ok(store)
    }

    fn persist_genomic_region_store(
        &mut self,
        mut store: gp::GenomicRegionStore,
    ) -> Result<(), EngineError> {
        store
            .sets
            .sort_by(|left, right| left.set_id.cmp(&right.set_id));
        for set in &store.sets {
            validate_set_digest(set)?;
        }
        let value = serde_json::to_value(store).map_err(|error| {
            region_error(
                ErrorCode::Internal,
                format!("could not encode genomic-region store: {error}"),
            )
        })?;
        self.state
            .metadata
            .insert(GENOMIC_REGION_SETS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    /// Return the portable saved-region store without mutating project state.
    pub fn genomic_region_store_snapshot(&self) -> Result<gp::GenomicRegionStore, EngineError> {
        self.genomic_region_store()
    }

    /// Format a saved ROI for human communication with an explicit coordinate convention.
    pub fn genomic_region_human_copy(region: &gp::GenomicRegionOfInterest) -> String {
        human_copy(region)
    }

    /// Format one BED6 row; BED coordinates are always 0-based and half-open.
    pub fn genomic_region_bed_row(region: &gp::GenomicRegionOfInterest) -> String {
        bed_row(region)
    }

    /// Serialize one complete canonical ROI for lossless clipboard interchange.
    pub fn genomic_region_canonical_json(
        region: &gp::GenomicRegionOfInterest,
    ) -> Result<String, EngineError> {
        serde_json::to_string_pretty(region).map_err(|error| {
            region_error(
                ErrorCode::Internal,
                format!("could not serialize canonical ROI JSON: {error}"),
            )
        })
    }

    pub(crate) fn project_genomic_region_sets_for_locus(
        &self,
        seq_id: &str,
        set_ids: &[String],
        locus_local_start_1based: usize,
        locus_local_end_1based: usize,
    ) -> Result<(Vec<gp::GeneLocusSavedRegionOverlayRow>, Vec<String>), EngineError> {
        if set_ids.is_empty() {
            return Ok((vec![], vec![]));
        }
        let store = self.genomic_region_store()?;
        let mut rows = vec![];
        let mut warnings = vec![];
        let mut requested = BTreeSet::new();
        for set_id in set_ids {
            if !requested.insert(set_id) {
                continue;
            }
            let Some(set) = store.sets.iter().find(|set| &set.set_id == set_id) else {
                warnings.push(format!(
                    "saved genomic region set '{set_id}' is unavailable"
                ));
                continue;
            };
            for region in &set.regions {
                let projection = match self.local_projection_for_interval(seq_id, &region.interval)
                {
                    Ok(projection) => projection,
                    Err(error) => {
                        warnings.push(format!(
                            "saved region '{}:{}' is not compatible with locus sequence '{}': {}",
                            set_id, region.region_id, seq_id, error.message
                        ));
                        continue;
                    }
                };
                let local_start_1based = projection.local_start_0based.saturating_add(1) as usize;
                let local_end_1based = projection.local_end_0based_exclusive as usize;
                if local_end_1based < locus_local_start_1based
                    || local_start_1based > locus_local_end_1based
                {
                    continue;
                }
                let availability =
                    if region.evidence.iter().any(|item| {
                        item.availability == gp::GenomicRegionEvidenceAvailability::Stale
                    }) {
                        gp::GenomicRegionEvidenceAvailability::Stale
                    } else if region.evidence.iter().any(|item| {
                        item.availability == gp::GenomicRegionEvidenceAvailability::Unavailable
                    }) {
                        gp::GenomicRegionEvidenceAvailability::Unavailable
                    } else if region.evidence.iter().any(|item| {
                        item.availability == gp::GenomicRegionEvidenceAvailability::Unverified
                    }) {
                        gp::GenomicRegionEvidenceAvailability::Unverified
                    } else {
                        gp::GenomicRegionEvidenceAvailability::Available
                    };
                rows.push(gp::GeneLocusSavedRegionOverlayRow {
                    set_id: set.set_id.clone(),
                    set_content_sha256: set.content_sha256.clone(),
                    region_id: region.region_id.clone(),
                    label: region.label.clone(),
                    purpose: region.purpose,
                    display_color_hex: region.display_color_hex.clone(),
                    selection_method: region.selection_method,
                    evidence_availability: availability,
                    local_start_1based,
                    local_end_1based,
                    genomic_start_0based: region.interval.start_0based,
                    genomic_end_0based_exclusive: region.interval.end_0based_exclusive,
                    genomic_strand: region.interval.strand,
                    region_content_sha256: region.content_sha256.clone(),
                    evidence_ids: region
                        .evidence
                        .iter()
                        .map(|item| item.evidence_id.clone())
                        .collect(),
                });
            }
        }
        rows.sort_by(|left, right| {
            (
                left.local_start_1based,
                left.local_end_1based,
                &left.set_id,
                &left.region_id,
            )
                .cmp(&(
                    right.local_start_1based,
                    right.local_end_1based,
                    &right.set_id,
                    &right.region_id,
                ))
        });
        Ok((rows, warnings))
    }

    fn anchor_reference(
        anchor: &SequenceGenomeAnchorSummary,
        reference_override: Option<gp::GenomicRegionReference>,
    ) -> Result<gp::GenomicRegionReference, EngineError> {
        let fallback = gp::GenomicRegionReference {
            assembly_name: anchor.genome_id.clone(),
            contig_name: anchor.chromosome.clone(),
            ..Default::default()
        };
        let reference = reference_override.unwrap_or(fallback);
        validate_reference(&reference)?;
        let anchor_assembly_matches = reference.assembly_name == anchor.genome_id
            || reference.assembly_accession.as_deref() == Some(anchor.genome_id.as_str());
        let contig_matches = reference.contig_name == anchor.chromosome
            || reference.contig_accession.as_deref() == Some(anchor.chromosome.as_str())
            || reference
                .contig_aliases
                .iter()
                .any(|alias| alias == &anchor.chromosome);
        if !anchor_assembly_matches || !contig_matches {
            return Err(region_error(
                ErrorCode::InvalidInput,
                format!(
                    "reference override {} / {} does not match exact sequence anchor {} / {}",
                    reference.assembly_name,
                    reference.contig_name,
                    anchor.genome_id,
                    anchor.chromosome
                ),
            ));
        }
        Ok(reference)
    }

    fn interval_and_projection_from_local(
        &self,
        seq_id: &str,
        local_start_0based: u64,
        local_end_0based_exclusive: u64,
        local_strand: gp::GenomicRegionStrand,
        reference_override: Option<gp::GenomicRegionReference>,
    ) -> Result<(gp::GenomicRegionInterval, gp::GenomicRegionLocalProjection), EngineError> {
        let dna = self.state.sequences.get(seq_id).ok_or_else(|| {
            region_error(
                ErrorCode::NotFound,
                format!("Sequence '{seq_id}' not found"),
            )
        })?;
        let sequence_len = u64::try_from(dna.len()).map_err(|_| {
            region_error(ErrorCode::InvalidInput, "sequence length does not fit u64")
        })?;
        if local_start_0based >= local_end_0based_exclusive
            || local_end_0based_exclusive > sequence_len
        {
            return Err(region_error(
                ErrorCode::InvalidInput,
                format!(
                    "local ROI {}..{} is outside sequence '{}' ({} bp)",
                    local_start_0based, local_end_0based_exclusive, seq_id, sequence_len
                ),
            ));
        }
        let anchor = self.sequence_genome_anchor_summary(seq_id)?;
        let reference = Self::anchor_reference(&anchor, reference_override)?;
        let anchor_start_0based = u64::try_from(anchor.start_1based.saturating_sub(1))
            .map_err(|_| region_error(ErrorCode::InvalidInput, "anchor start does not fit u64"))?;
        let anchor_end_0based_exclusive = u64::try_from(anchor.end_1based)
            .map_err(|_| region_error(ErrorCode::InvalidInput, "anchor end does not fit u64"))?;
        let anchor_strand = strand_from_char(anchor.strand);
        let (start_0based, end_0based_exclusive, genomic_strand) =
            if anchor_strand == gp::GenomicRegionStrand::Minus {
                (
                    anchor_end_0based_exclusive
                        .checked_sub(local_end_0based_exclusive)
                        .ok_or_else(|| {
                            region_error(ErrorCode::InvalidInput, "local ROI underflows anchor")
                        })?,
                    anchor_end_0based_exclusive
                        .checked_sub(local_start_0based)
                        .ok_or_else(|| {
                            region_error(ErrorCode::InvalidInput, "local ROI underflows anchor")
                        })?,
                    flip_strand(local_strand),
                )
            } else {
                (
                    anchor_start_0based
                        .checked_add(local_start_0based)
                        .ok_or_else(|| {
                            region_error(ErrorCode::InvalidInput, "genomic ROI start overflows")
                        })?,
                    anchor_start_0based
                        .checked_add(local_end_0based_exclusive)
                        .ok_or_else(|| {
                            region_error(ErrorCode::InvalidInput, "genomic ROI end overflows")
                        })?,
                    local_strand,
                )
            };
        if start_0based < anchor_start_0based || end_0based_exclusive > anchor_end_0based_exclusive
        {
            return Err(region_error(
                ErrorCode::InvalidInput,
                "projected genomic ROI falls outside the exact sequence anchor",
            ));
        }
        Ok((
            gp::GenomicRegionInterval {
                reference,
                start_0based,
                end_0based_exclusive,
                strand: genomic_strand,
                coordinate_convention: gp::GenomicRegionCoordinateConvention::ZeroBasedHalfOpen,
            },
            gp::GenomicRegionLocalProjection {
                seq_id: seq_id.to_string(),
                sequence_sha256: sha256_prefixed_bytes(dna.forward_bytes()),
                source_genome_id: anchor.genome_id,
                anchor_start_1based: anchor.start_1based as u64,
                anchor_end_1based: anchor.end_1based as u64,
                anchor_strand,
                local_start_0based,
                local_end_0based_exclusive,
                local_strand,
                status: gp::GenomicRegionLocalProjectionStatus::Current,
            },
        ))
    }

    fn local_projection_for_interval(
        &self,
        seq_id: &str,
        interval: &gp::GenomicRegionInterval,
    ) -> Result<gp::GenomicRegionLocalProjection, EngineError> {
        let anchor = self.sequence_genome_anchor_summary(seq_id)?;
        let reference = Self::anchor_reference(&anchor, Some(interval.reference.clone()))?;
        if !references_match(&reference, &interval.reference) {
            return Err(region_error(
                ErrorCode::InvalidInput,
                "ROI reference does not match the requested sequence anchor",
            ));
        }
        let anchor_start_0based = anchor.start_1based.saturating_sub(1) as u64;
        let anchor_end_0based_exclusive = anchor.end_1based as u64;
        if interval.start_0based < anchor_start_0based
            || interval.end_0based_exclusive > anchor_end_0based_exclusive
        {
            return Err(region_error(
                ErrorCode::InvalidInput,
                "ROI lies outside the requested sequence anchor",
            ));
        }
        let anchor_strand = strand_from_char(anchor.strand);
        let (local_start, local_end, local_strand) =
            if anchor_strand == gp::GenomicRegionStrand::Minus {
                (
                    anchor_end_0based_exclusive - interval.end_0based_exclusive,
                    anchor_end_0based_exclusive - interval.start_0based,
                    flip_strand(interval.strand),
                )
            } else {
                (
                    interval.start_0based - anchor_start_0based,
                    interval.end_0based_exclusive - anchor_start_0based,
                    interval.strand,
                )
            };
        let dna = self.state.sequences.get(seq_id).ok_or_else(|| {
            region_error(
                ErrorCode::NotFound,
                format!("Sequence '{seq_id}' not found"),
            )
        })?;
        Ok(gp::GenomicRegionLocalProjection {
            seq_id: seq_id.to_string(),
            sequence_sha256: sha256_prefixed_bytes(dna.forward_bytes()),
            source_genome_id: anchor.genome_id,
            anchor_start_1based: anchor.start_1based as u64,
            anchor_end_1based: anchor.end_1based as u64,
            anchor_strand,
            local_start_0based: local_start,
            local_end_0based_exclusive: local_end,
            local_strand,
            status: gp::GenomicRegionLocalProjectionStatus::Current,
        })
    }

    fn refresh_projection_status(
        &self,
        region: &mut gp::GenomicRegionOfInterest,
    ) -> Result<(), EngineError> {
        let Some(projection) = region.local_projection.as_mut() else {
            return Ok(());
        };
        let Some(dna) = self.state.sequences.get(&projection.seq_id) else {
            projection.status = gp::GenomicRegionLocalProjectionStatus::SequenceUnavailable;
            return Ok(());
        };
        if sha256_prefixed_bytes(dna.forward_bytes()) != projection.sequence_sha256 {
            projection.status = gp::GenomicRegionLocalProjectionStatus::SequenceDigestMismatch;
            return Ok(());
        }
        let Ok(anchor) = self.sequence_genome_anchor_summary(&projection.seq_id) else {
            projection.status = gp::GenomicRegionLocalProjectionStatus::AnchorMismatch;
            return Ok(());
        };
        if anchor.genome_id != projection.source_genome_id
            || anchor.start_1based as u64 != projection.anchor_start_1based
            || anchor.end_1based as u64 != projection.anchor_end_1based
            || strand_from_char(anchor.strand) != projection.anchor_strand
        {
            projection.status = gp::GenomicRegionLocalProjectionStatus::AnchorMismatch;
            return Ok(());
        }
        projection.status = gp::GenomicRegionLocalProjectionStatus::Current;
        Ok(())
    }

    fn insert_region(
        &mut self,
        set_id: &str,
        set_label: Option<String>,
        region: gp::GenomicRegionOfInterest,
        collision: gp::GenomicRegionCollisionPolicy,
    ) -> Result<gp::GenomicRegionSet, EngineError> {
        valid_identifier(set_id, "set_id")?;
        let mut store = self.genomic_region_store()?;
        let set_index = store.sets.iter().position(|set| set.set_id == set_id);
        if set_index.is_none() {
            store.sets.push(gp::GenomicRegionSet {
                schema: gp::GENOMIC_REGION_SET_SCHEMA.to_string(),
                set_id: set_id.to_string(),
                label: set_label.clone(),
                description: None,
                regions: vec![],
                content_sha256: String::new(),
            });
        }
        let index = store
            .sets
            .iter()
            .position(|set| set.set_id == set_id)
            .ok_or_else(|| region_error(ErrorCode::Internal, "new region set was not retained"))?;
        let set = &mut store.sets[index];
        if set.label.is_none() {
            set.label = set_label;
        }
        if let Some(existing) = set
            .regions
            .iter()
            .position(|item| item.region_id == region.region_id)
        {
            if collision == gp::GenomicRegionCollisionPolicy::Reject {
                return Err(region_error(
                    ErrorCode::InvalidInput,
                    format!(
                        "region '{}' already exists in set '{}'; request replace explicitly",
                        region.region_id, set_id
                    ),
                ));
            }
            set.regions[existing] = region;
        } else {
            set.regions.push(region);
        }
        recompute_set_digest(set)?;
        validate_set_digest(set)?;
        let saved = set.clone();
        self.persist_genomic_region_store(store)?;
        Ok(saved)
    }

    fn create_region(
        &mut self,
        mut request: gp::GenomicRegionCreateRequest,
    ) -> Result<gp::GenomicRegionOperationReport, EngineError> {
        validate_interval(&request.interval)?;
        request.display_color_hex = normalize_display_color(request.display_color_hex)?;
        if let Some(supplied) = request.local_projection.as_ref() {
            let expected =
                self.local_projection_for_interval(&supplied.seq_id, &request.interval)?;
            let mut expected_identity = expected.clone();
            expected_identity.status = supplied.status;
            if &expected_identity != supplied {
                return Err(region_error(
                    ErrorCode::InvalidInput,
                    "supplied local projection does not match the current sequence digest, genome anchor, interval, strand, or local coordinates",
                ));
            }
        }
        if request.selection_method == gp::GenomicRegionSelectionMethod::Derived {
            return Err(region_error(
                ErrorCode::InvalidInput,
                "use DeriveGenomicRegion for derived geometry",
            ));
        }
        let mut region = gp::GenomicRegionOfInterest {
            schema: gp::GENOMIC_REGION_OF_INTEREST_SCHEMA.to_string(),
            region_id: request.region_id.unwrap_or_default(),
            label: request.label,
            description: request.description,
            interval: request.interval,
            local_projection: request.local_projection,
            purpose: request.purpose,
            display_color_hex: request.display_color_hex,
            selection_method: request.selection_method,
            evidence: request.evidence,
            derivation: None,
            created_at_unix_ms: request.created_at_unix_ms,
            notes: request.notes,
            identity_sha256: String::new(),
            content_sha256: String::new(),
        };
        recompute_region_digests(&mut region)?;
        valid_identifier(&region.region_id, "region_id")?;
        let set = self.insert_region(
            &request.set_id,
            request.set_label,
            region.clone(),
            request.collision_policy,
        )?;
        report_for_region("create", Some(set), Some(region))
    }

    fn capture_region(
        &mut self,
        request: gp::GenomicRegionCaptureRequest,
    ) -> Result<gp::GenomicRegionOperationReport, EngineError> {
        let (interval, local_projection, selection_method, evidence) = match request.source {
            gp::GenomicRegionCaptureSource::SequenceSelection {
                seq_id,
                local_start_0based,
                local_end_0based_exclusive,
                strand,
                reference_override,
            } => {
                let (interval, projection) = self.interval_and_projection_from_local(
                    &seq_id,
                    local_start_0based,
                    local_end_0based_exclusive,
                    strand,
                    reference_override,
                )?;
                (
                    interval,
                    Some(projection),
                    gp::GenomicRegionSelectionMethod::ManualSpan,
                    vec![],
                )
            }
            gp::GenomicRegionCaptureSource::SequenceFeature {
                seq_id,
                feature_id,
                reference_override,
            } => {
                let dna = self.state.sequences.get(&seq_id).ok_or_else(|| {
                    region_error(
                        ErrorCode::NotFound,
                        format!("Sequence '{seq_id}' not found"),
                    )
                })?;
                let feature = dna.features().get(feature_id).ok_or_else(|| {
                    region_error(
                        ErrorCode::NotFound,
                        format!("Feature {feature_id} not found on '{seq_id}'"),
                    )
                })?;
                let (start, end) = feature.location.find_bounds().map_err(|error| {
                    region_error(
                        ErrorCode::InvalidInput,
                        format!("Feature {feature_id} has no bounded interval: {error}"),
                    )
                })?;
                let start = u64::try_from(start).map_err(|_| {
                    region_error(ErrorCode::InvalidInput, "feature start is negative")
                })?;
                let end = u64::try_from(end).map_err(|_| {
                    region_error(ErrorCode::InvalidInput, "feature end is negative")
                })?;
                let strand = if feature_is_reverse(feature) {
                    gp::GenomicRegionStrand::Minus
                } else {
                    gp::GenomicRegionStrand::Plus
                };
                let (interval, projection) = self.interval_and_projection_from_local(
                    &seq_id,
                    start,
                    end,
                    strand,
                    reference_override,
                )?;
                let sequence_sha = projection.sequence_sha256.clone();
                let evidence = gp::GenomicRegionEvidenceReference {
                    evidence_id: format!("feature:{seq_id}:{feature_id}"),
                    source_kind: "sequence_feature".to_string(),
                    source_id: seq_id,
                    source_sha256: Some(sequence_sha),
                    feature_or_window_id: Some(feature_id.to_string()),
                    feature_type: Some(feature.kind.to_string()),
                    availability: gp::GenomicRegionEvidenceAvailability::Available,
                    evidence_statement:
                        "Existing sequence annotation supplied the selected interval geometry."
                            .to_string(),
                    non_claims: vec![
                        "Annotation does not by itself establish regulatory activity or causality."
                            .to_string(),
                    ],
                    ..Default::default()
                };
                (
                    interval,
                    Some(projection),
                    gp::GenomicRegionSelectionMethod::ExistingInterval,
                    vec![evidence],
                )
            }
            gp::GenomicRegionCaptureSource::CutrunSupportWindow {
                report,
                window_id,
                reference,
            } => {
                validate_reference(&reference)?;
                let window = report
                    .support_windows
                    .iter()
                    .find(|window| window.window_id == window_id)
                    .ok_or_else(|| {
                        region_error(
                            ErrorCode::NotFound,
                            format!("CUT&RUN support window '{window_id}' not found"),
                        )
                    })?;
                if report.seq_id.trim().is_empty() {
                    return Err(region_error(
                        ErrorCode::InvalidInput,
                        "CUT&RUN support report has no sequence id",
                    ));
                }
                let interval = gp::GenomicRegionInterval {
                    reference,
                    start_0based: window.genomic_start_1based.saturating_sub(1) as u64,
                    end_0based_exclusive: window.genomic_end_1based as u64,
                    strand: gp::GenomicRegionStrand::Unstranded,
                    coordinate_convention: gp::GenomicRegionCoordinateConvention::ZeroBasedHalfOpen,
                };
                validate_interval(&interval)?;
                let projection =
                    Some(self.local_projection_for_interval(&report.seq_id, &interval)?);
                let report_sha = sha256_prefixed_str(&serialize_for_digest(
                    report.as_ref(),
                    "CUT&RUN support report",
                )?);
                let mut evidence = report
                    .evidence_sources
                    .iter()
                    .enumerate()
                    .map(|(index, source)| {
                        let source_kind = match source.source_kind {
                            gp::CutRunRegulatoryEvidenceSourceKind::Dataset => "cutrun_dataset",
                            gp::CutRunRegulatoryEvidenceSourceKind::ReadReport => {
                                "cutrun_read_report"
                            }
                        };
                        gp::GenomicRegionEvidenceReference {
                            evidence_id: format!(
                                "cutrun:{}:{window_id}:{source_kind}:{}",
                                report.seq_id,
                                index + 1
                            ),
                            source_kind: source_kind.to_string(),
                            source_id: source.source_id.clone(),
                            source_sha256: Some(report_sha.clone()),
                            report_id: source.report_id.clone(),
                            feature_or_window_id: Some(window_id.clone()),
                            assay: Some("CUT&RUN".to_string()),
                            factor: source.target_factor.clone(),
                            support_strength: Some(
                                format!("{:?}", window.support_strength).to_ascii_lowercase(),
                            ),
                            max_signal_value: window.max_signal_value,
                            availability: gp::GenomicRegionEvidenceAvailability::Available,
                            evidence_statement: "CUT&RUN enrichment/support selected this genomic window as occupancy evidence.".to_string(),
                            non_claims: vec![
                                "CUT&RUN enrichment is not biochemical affinity.".to_string(),
                                "Occupancy does not by itself establish causal regulatory activity.".to_string(),
                            ],
                            ..Default::default()
                        }
                    })
                    .collect::<Vec<_>>();
                if evidence.is_empty() {
                    evidence.push(gp::GenomicRegionEvidenceReference {
                        evidence_id: format!("cutrun:{}:{window_id}:support_report", report.seq_id),
                        source_kind: "cutrun_support_report".to_string(),
                        source_id: report.seq_id.clone(),
                        source_sha256: Some(report_sha),
                        feature_or_window_id: Some(window_id),
                        assay: Some("CUT&RUN".to_string()),
                        support_strength: Some(
                            format!("{:?}", window.support_strength).to_ascii_lowercase(),
                        ),
                        max_signal_value: window.max_signal_value,
                        availability: gp::GenomicRegionEvidenceAvailability::Available,
                        evidence_statement: "CUT&RUN enrichment/support selected this genomic window as occupancy evidence.".to_string(),
                        non_claims: vec![
                            "CUT&RUN enrichment is not biochemical affinity.".to_string(),
                            "Occupancy does not by itself establish causal regulatory activity."
                                .to_string(),
                        ],
                        ..Default::default()
                    });
                }
                (
                    interval,
                    projection,
                    gp::GenomicRegionSelectionMethod::CutrunSupportWindow,
                    evidence,
                )
            }
            gp::GenomicRegionCaptureSource::EnsemblRegulatoryFeature {
                report,
                feature_id,
                interval_kind,
            } => {
                if !report.content_identity_verified {
                    return Err(region_error(
                        ErrorCode::InvalidInput,
                        "Ensembl Regulation report content identity is not verified",
                    ));
                }
                let row = report
                    .rows
                    .iter()
                    .find(|row| row.interval.feature_id == feature_id)
                    .ok_or_else(|| {
                        region_error(
                            ErrorCode::NotFound,
                            format!("Ensembl regulatory feature '{feature_id}' not found"),
                        )
                    })?;
                let (start, end) = match interval_kind {
                    gp::GenomicRegionEnsemblIntervalKind::Core => {
                        (row.interval.start_0based, row.interval.end_0based_exclusive)
                    }
                    gp::GenomicRegionEnsemblIntervalKind::Extended => (
                        row.interval.extended_start_0based.ok_or_else(|| {
                            region_error(
                                ErrorCode::InvalidInput,
                                format!("feature '{feature_id}' has no extended interval"),
                            )
                        })?,
                        row.interval.extended_end_0based_exclusive.ok_or_else(|| {
                            region_error(
                                ErrorCode::InvalidInput,
                                format!("feature '{feature_id}' has no extended interval"),
                            )
                        })?,
                    ),
                };
                let interval = gp::GenomicRegionInterval {
                    reference: gp::GenomicRegionReference {
                        species_scientific_name: Some(
                            report.source.species_scientific_name.clone(),
                        ),
                        taxon_id: Some(report.source.taxon_id),
                        assembly_name: report.source.assembly_name.clone(),
                        assembly_accession: Some(report.source.assembly_accession.clone())
                            .filter(|value| !value.trim().is_empty()),
                        contig_name: row.interval.chromosome.clone(),
                        contig_accession: None,
                        contig_aliases: vec![],
                    },
                    start_0based: start,
                    end_0based_exclusive: end,
                    strand: strand_from_char(row.interval.strand),
                    coordinate_convention: gp::GenomicRegionCoordinateConvention::ZeroBasedHalfOpen,
                };
                validate_interval(&interval)?;
                let projection =
                    Some(self.local_projection_for_interval(&report.seq_id, &interval)?);
                let canonical_url = if report
                    .source
                    .feature_page_url_template
                    .matches("{FEATURE_ID}")
                    .count()
                    == 1
                {
                    Some(
                        report
                            .source
                            .feature_page_url_template
                            .replace("{FEATURE_ID}", &feature_id),
                    )
                } else {
                    None
                };
                let report_sha = sha256_prefixed_str(&serialize_for_digest(
                    report.as_ref(),
                    "Ensembl Regulation overlap report",
                )?);
                let evidence = gp::GenomicRegionEvidenceReference {
                    evidence_id: format!("ensembl_regulation:{feature_id}:{interval_kind:?}"),
                    source_kind: "ensembl_regulatory_feature".to_string(),
                    source_id: report.source.source_id.clone(),
                    source_release: Some(report.source.annotation_release.clone()),
                    source_sha256: Some(report_sha),
                    report_id: Some(report.report_id.clone()),
                    feature_or_window_id: Some(feature_id),
                    feature_type: Some(row.interval.feature_type.clone()),
                    provider_url: canonical_url,
                    associated_gene_ids: row.interval.associated_gene_ids.clone(),
                    associated_gene_names: row.interval.associated_gene_names.clone(),
                    availability: gp::GenomicRegionEvidenceAvailability::Available,
                    evidence_statement: report.evidence_statement.clone(),
                    non_claims: report.non_claims.clone(),
                    ..Default::default()
                };
                (
                    interval,
                    projection,
                    gp::GenomicRegionSelectionMethod::EnsemblRegulatoryFeature,
                    vec![evidence],
                )
            }
            gp::GenomicRegionCaptureSource::GeneLocusEnsemblRegulatoryFeature {
                row,
                source_binding,
                reference,
                seq_id,
                interval_kind,
                evidence_statement,
                non_claims,
            } => {
                if !source_binding.content_identity_verified {
                    return Err(region_error(
                        ErrorCode::InvalidInput,
                        "gene-locus Ensembl Regulation source identity is not verified",
                    ));
                }
                if row.source_id != source_binding.source.source_id
                    || row.annotation_release != source_binding.source.annotation_release
                    || row.assembly_name != source_binding.source.assembly_name
                {
                    return Err(region_error(
                        ErrorCode::InvalidInput,
                        "gene-locus Ensembl row does not match its source binding",
                    ));
                }
                validate_reference(&reference)?;
                if reference.assembly_name != row.assembly_name
                    || reference
                        .assembly_accession
                        .as_deref()
                        .is_some_and(|value| value != row.assembly_accession.as_str())
                {
                    return Err(region_error(
                        ErrorCode::InvalidInput,
                        "gene-locus Ensembl row does not match the explicit ROI reference",
                    ));
                }
                let (start_1based, end_1based) = match interval_kind {
                    gp::GenomicRegionEnsemblIntervalKind::Core => {
                        (row.core_genomic_start_1based, row.core_genomic_end_1based)
                    }
                    gp::GenomicRegionEnsemblIntervalKind::Extended => (
                        row.extended_genomic_start_1based.ok_or_else(|| {
                            region_error(
                                ErrorCode::InvalidInput,
                                format!("feature '{}' has no extended interval", row.feature_id),
                            )
                        })?,
                        row.extended_genomic_end_1based.ok_or_else(|| {
                            region_error(
                                ErrorCode::InvalidInput,
                                format!("feature '{}' has no extended interval", row.feature_id),
                            )
                        })?,
                    ),
                };
                let interval = gp::GenomicRegionInterval {
                    reference,
                    start_0based: u64::try_from(start_1based.saturating_sub(1)).map_err(|_| {
                        region_error(
                            ErrorCode::InvalidInput,
                            "Ensembl interval start does not fit u64",
                        )
                    })?,
                    end_0based_exclusive: u64::try_from(end_1based).map_err(|_| {
                        region_error(
                            ErrorCode::InvalidInput,
                            "Ensembl interval end does not fit u64",
                        )
                    })?,
                    strand: strand_from_char(row.genomic_strand.chars().next()),
                    coordinate_convention: gp::GenomicRegionCoordinateConvention::ZeroBasedHalfOpen,
                };
                validate_interval(&interval)?;
                let projection = seq_id
                    .as_deref()
                    .map(|seq_id| self.local_projection_for_interval(seq_id, &interval))
                    .transpose()?;
                let binding_sha = sha256_prefixed_str(&serialize_for_digest(
                    &source_binding,
                    "gene-locus Ensembl source binding",
                )?);
                let evidence = gp::GenomicRegionEvidenceReference {
                    evidence_id: format!("ensembl_regulation:{}:{interval_kind:?}", row.feature_id),
                    source_kind: "ensembl_regulatory_feature".to_string(),
                    source_id: row.source_id.clone(),
                    source_release: Some(row.annotation_release.clone()),
                    source_sha256: Some(binding_sha),
                    report_id: Some(source_binding.overlap_report_id),
                    feature_or_window_id: Some(row.feature_id),
                    feature_type: Some(row.feature_type),
                    provider_url: Some(row.canonical_feature_url)
                        .filter(|value| !value.trim().is_empty()),
                    associated_gene_ids: row.associated_gene_ids,
                    associated_gene_names: row.associated_gene_names,
                    availability: gp::GenomicRegionEvidenceAvailability::Available,
                    evidence_statement,
                    non_claims,
                    ..Default::default()
                };
                (
                    interval,
                    projection,
                    gp::GenomicRegionSelectionMethod::EnsemblRegulatoryFeature,
                    vec![evidence],
                )
            }
            gp::GenomicRegionCaptureSource::ProviderAnnotation { interval, evidence } => {
                validate_interval(&interval)?;
                (
                    interval,
                    None,
                    gp::GenomicRegionSelectionMethod::ProviderAnnotation,
                    vec![evidence],
                )
            }
        };
        self.create_region(gp::GenomicRegionCreateRequest {
            set_id: request.set_id,
            set_label: request.set_label,
            region_id: request.region_id,
            label: request.label,
            description: request.description,
            interval,
            local_projection,
            purpose: request.purpose,
            display_color_hex: request.display_color_hex,
            selection_method,
            evidence,
            created_at_unix_ms: request.created_at_unix_ms,
            notes: request.notes,
            collision_policy: request.collision_policy,
        })
    }

    fn list_regions(
        &self,
        request: gp::GenomicRegionListRequest,
    ) -> Result<gp::GenomicRegionOperationReport, EngineError> {
        let store = self.genomic_region_store()?;
        let selected_set = request
            .set_id
            .as_ref()
            .and_then(|set_id| store.sets.iter().find(|set| &set.set_id == set_id))
            .cloned();
        let sets = store
            .sets
            .iter()
            .filter(|set| request.set_id.as_ref().is_none_or(|id| id == &set.set_id))
            .map(|set| gp::GenomicRegionSetSummary {
                set_id: set.set_id.clone(),
                label: set.label.clone(),
                region_count: set.regions.len(),
                content_sha256: set.content_sha256.clone(),
            })
            .collect::<Vec<_>>();
        if request.set_id.is_some() && sets.is_empty() {
            return Err(region_error(
                ErrorCode::NotFound,
                "genomic region set not found",
            ));
        }
        Ok(gp::GenomicRegionOperationReport {
            schema: gp::GENOMIC_REGION_OPERATION_REPORT_SCHEMA.to_string(),
            action: "list".to_string(),
            set: selected_set,
            sets,
            ..Default::default()
        })
    }

    fn inspect_region(
        &self,
        request: gp::GenomicRegionInspectRequest,
    ) -> Result<gp::GenomicRegionOperationReport, EngineError> {
        let store = self.genomic_region_store()?;
        let set = store
            .sets
            .iter()
            .find(|set| set.set_id == request.set_id)
            .ok_or_else(|| {
                region_error(
                    ErrorCode::NotFound,
                    format!("genomic region set '{}' not found", request.set_id),
                )
            })?;
        let mut region = set
            .regions
            .iter()
            .find(|region| region.region_id == request.region_id)
            .cloned()
            .ok_or_else(|| {
                region_error(
                    ErrorCode::NotFound,
                    format!("genomic region '{}' not found", request.region_id),
                )
            })?;
        self.refresh_projection_status(&mut region)?;
        recompute_region_digests(&mut region)?;
        let mut inspected_set = set.clone();
        if let Some(saved) = inspected_set
            .regions
            .iter_mut()
            .find(|saved| saved.region_id == request.region_id)
        {
            *saved = region.clone();
        }
        recompute_set_digest(&mut inspected_set)?;
        let mut report = report_for_region("inspect", Some(inspected_set), Some(region))?;
        if report
            .region
            .as_ref()
            .and_then(|region| region.local_projection.as_ref())
            .is_some_and(|projection| {
                projection.status != gp::GenomicRegionLocalProjectionStatus::Current
            })
        {
            report.warnings.push(
                "The saved canonical genomic interval remains available, but its local projection is not current."
                    .to_string(),
            );
        }
        Ok(report)
    }

    fn update_region(
        &mut self,
        request: gp::GenomicRegionUpdateRequest,
    ) -> Result<gp::GenomicRegionOperationReport, EngineError> {
        let mut store = self.genomic_region_store()?;
        let set = store
            .sets
            .iter_mut()
            .find(|set| set.set_id == request.set_id)
            .ok_or_else(|| region_error(ErrorCode::NotFound, "genomic region set not found"))?;
        let region = set
            .regions
            .iter_mut()
            .find(|region| region.region_id == request.region_id)
            .ok_or_else(|| region_error(ErrorCode::NotFound, "genomic region not found"))?;
        let identity_before = region.identity_sha256.clone();
        region.label = request.label;
        region.description = request.description;
        if request.display_color_hex.is_some() {
            region.display_color_hex = normalize_display_color(request.display_color_hex)?;
        }
        region.notes = request.notes;
        recompute_region_digests(region)?;
        if region.identity_sha256 != identity_before {
            return Err(region_error(
                ErrorCode::Internal,
                "presentation-only update changed genomic-region identity",
            ));
        }
        let updated = region.clone();
        recompute_set_digest(set)?;
        let saved_set = set.clone();
        self.persist_genomic_region_store(store)?;
        report_for_region("update_presentation", Some(saved_set), Some(updated))
    }

    fn derive_region(
        &mut self,
        request: gp::GenomicRegionDeriveRequest,
    ) -> Result<gp::GenomicRegionOperationReport, EngineError> {
        if request.parent_region_ids.len() < 2 {
            return Err(region_error(
                ErrorCode::InvalidInput,
                "derived regions require at least two parent IDs",
            ));
        }
        let store = self.genomic_region_store()?;
        let set = store
            .sets
            .iter()
            .find(|set| set.set_id == request.set_id)
            .ok_or_else(|| region_error(ErrorCode::NotFound, "genomic region set not found"))?;
        let mut parent_ids = BTreeSet::new();
        let mut parents = vec![];
        for id in &request.parent_region_ids {
            if !parent_ids.insert(id) {
                return Err(region_error(
                    ErrorCode::InvalidInput,
                    format!("duplicate parent region id '{id}'"),
                ));
            }
            parents.push(
                set.regions
                    .iter()
                    .find(|region| &region.region_id == id)
                    .cloned()
                    .ok_or_else(|| {
                        region_error(
                            ErrorCode::NotFound,
                            format!("parent region '{id}' not found"),
                        )
                    })?,
            );
        }
        let first = &parents[0];
        if parents
            .iter()
            .skip(1)
            .any(|parent| !references_match(&first.interval.reference, &parent.interval.reference))
        {
            return Err(region_error(
                ErrorCode::InvalidInput,
                "parent regions have incompatible assembly or contig identities; no liftover is implicit",
            ));
        }
        let start = match request.method {
            gp::GenomicRegionDerivationMethod::Intersection => parents
                .iter()
                .map(|parent| parent.interval.start_0based)
                .max()
                .unwrap_or_default(),
            gp::GenomicRegionDerivationMethod::Union | gp::GenomicRegionDerivationMethod::Hull => {
                parents
                    .iter()
                    .map(|parent| parent.interval.start_0based)
                    .min()
                    .unwrap_or_default()
            }
        };
        let end = match request.method {
            gp::GenomicRegionDerivationMethod::Intersection => parents
                .iter()
                .map(|parent| parent.interval.end_0based_exclusive)
                .min()
                .unwrap_or_default(),
            gp::GenomicRegionDerivationMethod::Union | gp::GenomicRegionDerivationMethod::Hull => {
                parents
                    .iter()
                    .map(|parent| parent.interval.end_0based_exclusive)
                    .max()
                    .unwrap_or_default()
            }
        };
        if request.method == gp::GenomicRegionDerivationMethod::Intersection && start >= end {
            return Err(region_error(
                ErrorCode::InvalidInput,
                "parent regions have no non-empty intersection",
            ));
        }
        if request.method == gp::GenomicRegionDerivationMethod::Union {
            let mut spans = parents
                .iter()
                .map(|parent| {
                    (
                        parent.interval.start_0based,
                        parent.interval.end_0based_exclusive,
                    )
                })
                .collect::<Vec<_>>();
            spans.sort_unstable();
            let mut covered_end = spans[0].1;
            for (next_start, next_end) in spans.into_iter().skip(1) {
                if next_start > covered_end {
                    return Err(region_error(
                        ErrorCode::InvalidInput,
                        "union parents are disjoint; use the explicit hull method to bridge gaps",
                    ));
                }
                covered_end = covered_end.max(next_end);
            }
        }
        let strand = if parents
            .iter()
            .all(|parent| parent.interval.strand == first.interval.strand)
        {
            first.interval.strand
        } else {
            gp::GenomicRegionStrand::Unstranded
        };
        let interval = gp::GenomicRegionInterval {
            reference: first.interval.reference.clone(),
            start_0based: start,
            end_0based_exclusive: end,
            strand,
            coordinate_convention: gp::GenomicRegionCoordinateConvention::ZeroBasedHalfOpen,
        };
        let parent_projections = parents
            .iter()
            .map(|parent| parent.local_projection.as_ref())
            .collect::<Vec<_>>();
        let projection_binding = parent_projections
            .first()
            .and_then(|projection| *projection);
        let projections_are_shared = projection_binding.is_some()
            && parent_projections.iter().all(|projection| {
                projection.is_some_and(|candidate| {
                    let Some(expected) = projection_binding else {
                        return false;
                    };
                    candidate.seq_id == expected.seq_id
                        && candidate.sequence_sha256 == expected.sequence_sha256
                        && candidate.source_genome_id == expected.source_genome_id
                        && candidate.anchor_start_1based == expected.anchor_start_1based
                        && candidate.anchor_end_1based == expected.anchor_end_1based
                        && candidate.anchor_strand == expected.anchor_strand
                })
            });
        let local_projection = match (projections_are_shared, projection_binding) {
            (true, Some(binding)) => {
                Some(self.local_projection_for_interval(&binding.seq_id, &interval)?)
            }
            _ => None,
        };
        let mut evidence_by_id = std::collections::BTreeMap::new();
        for evidence in parents
            .iter()
            .flat_map(|parent| parent.evidence.iter().cloned())
        {
            if let Some(existing) = evidence_by_id.get(&evidence.evidence_id) {
                if existing != &evidence {
                    return Err(region_error(
                        ErrorCode::InvalidInput,
                        format!(
                            "parent regions carry conflicting records for evidence_id '{}'",
                            evidence.evidence_id
                        ),
                    ));
                }
            } else {
                evidence_by_id.insert(evidence.evidence_id.clone(), evidence);
            }
        }
        let evidence = evidence_by_id.into_values().collect::<Vec<_>>();
        let derivation = gp::GenomicRegionDerivation {
            method: request.method,
            parents: parents
                .iter()
                .map(|parent| gp::GenomicRegionParentReference {
                    region_id: parent.region_id.clone(),
                    identity_sha256: parent.identity_sha256.clone(),
                    content_sha256: parent.content_sha256.clone(),
                })
                .collect(),
            rule: match request.method {
                gp::GenomicRegionDerivationMethod::Union => {
                    "Contiguous/overlapping parent intervals were unioned without widening across a gap."
                }
                gp::GenomicRegionDerivationMethod::Intersection => {
                    "The non-empty geometric intersection of every parent interval was retained."
                }
                gp::GenomicRegionDerivationMethod::Hull => {
                    "The explicit minimum-to-maximum hull may include bases not present in any parent."
                }
            }
            .to_string(),
        };
        let mut region = gp::GenomicRegionOfInterest {
            schema: gp::GENOMIC_REGION_OF_INTEREST_SCHEMA.to_string(),
            region_id: request.region_id.unwrap_or_default(),
            label: request.label,
            description: request.description,
            interval,
            local_projection,
            purpose: request.purpose,
            display_color_hex: normalize_display_color(request.display_color_hex)?,
            selection_method: gp::GenomicRegionSelectionMethod::Derived,
            evidence,
            derivation: Some(derivation),
            created_at_unix_ms: request.created_at_unix_ms,
            notes: request.notes,
            identity_sha256: String::new(),
            content_sha256: String::new(),
        };
        recompute_region_digests(&mut region)?;
        let set = self.insert_region(
            &request.set_id,
            None,
            region.clone(),
            request.collision_policy,
        )?;
        let mut report = report_for_region("derive", Some(set), Some(region))?;
        if parent_projections
            .iter()
            .any(|projection| projection.is_some())
            && !projections_are_shared
        {
            report.warnings.push(
                "The parents do not share one exact sequence/anchor binding, so no local projection was inferred for the derived genomic interval."
                    .to_string(),
            );
        }
        Ok(report)
    }

    fn read_bounded(path: &str, max_bytes: u64) -> Result<Vec<u8>, EngineError> {
        let metadata = std::fs::metadata(path).map_err(|error| {
            region_error(
                ErrorCode::Io,
                format!("could not inspect '{path}': {error}"),
            )
        })?;
        if metadata.len() > max_bytes {
            return Err(region_error(
                ErrorCode::InvalidInput,
                format!("'{path}' exceeds the declared {max_bytes}-byte import limit"),
            ));
        }
        std::fs::read(path).map_err(|error| {
            region_error(ErrorCode::Io, format!("could not read '{path}': {error}"))
        })
    }

    fn parse_bare_bed(
        bytes: &[u8],
        reference: gp::GenomicRegionReference,
        set_id: String,
        max_rows: usize,
    ) -> Result<gp::GenomicRegionSet, EngineError> {
        validate_reference(&reference)?;
        valid_identifier(&set_id, "set_id_override")?;
        let text = std::str::from_utf8(bytes).map_err(|error| {
            region_error(
                ErrorCode::InvalidInput,
                format!("BED is not UTF-8: {error}"),
            )
        })?;
        let mut regions = vec![];
        let mut ids = BTreeSet::new();
        for (line_index, line) in text.lines().enumerate() {
            let trimmed = line.trim();
            if trimmed.is_empty()
                || trimmed.starts_with('#')
                || trimmed.starts_with("track ")
                || trimmed.starts_with("browser ")
            {
                continue;
            }
            if regions.len() >= max_rows {
                return Err(region_error(
                    ErrorCode::InvalidInput,
                    format!("BED exceeds the declared {max_rows}-row import limit"),
                ));
            }
            let columns = line.split('\t').collect::<Vec<_>>();
            if columns.len() < 3 {
                return Err(region_error(
                    ErrorCode::InvalidInput,
                    format!("BED line {} has fewer than three columns", line_index + 1),
                ));
            }
            let contig_matches = columns[0] == reference.contig_name
                || reference.contig_accession.as_deref() == Some(columns[0])
                || reference
                    .contig_aliases
                    .iter()
                    .any(|alias| alias == columns[0]);
            if !contig_matches {
                return Err(region_error(
                    ErrorCode::InvalidInput,
                    format!(
                        "BED line {} contig '{}' does not match explicit reference '{}'",
                        line_index + 1,
                        columns[0],
                        reference.contig_name
                    ),
                ));
            }
            let start = columns[1].parse::<u64>().map_err(|error| {
                region_error(
                    ErrorCode::InvalidInput,
                    format!("invalid BED start on line {}: {error}", line_index + 1),
                )
            })?;
            let end = columns[2].parse::<u64>().map_err(|error| {
                region_error(
                    ErrorCode::InvalidInput,
                    format!("invalid BED end on line {}: {error}", line_index + 1),
                )
            })?;
            let label = columns
                .get(3)
                .map(|value| value.to_string())
                .filter(|value| !value.trim().is_empty() && value != ".");
            let strand = match columns.get(5).copied() {
                Some("+") => gp::GenomicRegionStrand::Plus,
                Some("-") => gp::GenomicRegionStrand::Minus,
                Some(".") | None => gp::GenomicRegionStrand::Unstranded,
                Some(other) => {
                    return Err(region_error(
                        ErrorCode::InvalidInput,
                        format!("invalid BED strand '{other}' on line {}", line_index + 1),
                    ));
                }
            };
            let region_id = label
                .as_deref()
                .filter(|value| valid_identifier(value, "BED name").is_ok())
                .map(str::to_string)
                .unwrap_or_else(|| format!("bed_row_{:06}", regions.len() + 1));
            if !ids.insert(region_id.clone()) {
                return Err(region_error(
                    ErrorCode::InvalidInput,
                    format!("duplicate BED region id '{region_id}'"),
                ));
            }
            let mut region = gp::GenomicRegionOfInterest {
                schema: gp::GENOMIC_REGION_OF_INTEREST_SCHEMA.to_string(),
                region_id,
                label,
                interval: gp::GenomicRegionInterval {
                    reference: reference.clone(),
                    start_0based: start,
                    end_0based_exclusive: end,
                    strand,
                    coordinate_convention: gp::GenomicRegionCoordinateConvention::ZeroBasedHalfOpen,
                },
                purpose: gp::GenomicRegionPurpose::Other,
                selection_method: gp::GenomicRegionSelectionMethod::Imported,
                evidence: vec![gp::GenomicRegionEvidenceReference {
                    evidence_id: format!("bed_line:{}", line_index + 1),
                    source_kind: "imported_bed".to_string(),
                    source_id: "caller_supplied_bed".to_string(),
                    availability: gp::GenomicRegionEvidenceAvailability::Unverified,
                    evidence_statement: "Caller-supplied BED row supplied interval geometry; no biological evidence was inferred.".to_string(),
                    ..Default::default()
                }],
                ..Default::default()
            };
            validate_interval(&region.interval)?;
            recompute_region_digests(&mut region)?;
            regions.push(region);
        }
        let mut set = gp::GenomicRegionSet {
            schema: gp::GENOMIC_REGION_SET_SCHEMA.to_string(),
            set_id,
            label: None,
            description: Some(
                "Imported from bare BED with an explicit caller-supplied genome reference."
                    .to_string(),
            ),
            regions,
            content_sha256: String::new(),
        };
        recompute_set_digest(&mut set)?;
        Ok(set)
    }

    fn import_region_set(
        &mut self,
        request: gp::GenomicRegionImportRequest,
    ) -> Result<gp::GenomicRegionOperationReport, EngineError> {
        if request.max_bytes == 0 || request.max_rows == 0 {
            return Err(region_error(
                ErrorCode::InvalidInput,
                "import limits must both be greater than zero",
            ));
        }
        let bytes = Self::read_bounded(&request.path, request.max_bytes)?;
        let mut set = match request.format {
            gp::GenomicRegionImportFormat::Json => {
                let set: gp::GenomicRegionSet =
                    serde_json::from_slice(&bytes).map_err(|error| {
                        region_error(
                            ErrorCode::InvalidInput,
                            format!("could not parse genomic region-set JSON: {error}"),
                        )
                    })?;
                if set.regions.len() > request.max_rows {
                    return Err(region_error(
                        ErrorCode::InvalidInput,
                        "region-set JSON exceeds the declared row limit",
                    ));
                }
                validate_set_digest(&set)?;
                set
            }
            gp::GenomicRegionImportFormat::Bed => {
                if let Some(manifest_path) = request.manifest_path.as_deref() {
                    let manifest_bytes = Self::read_bounded(manifest_path, request.max_bytes)?;
                    let manifest: gp::GenomicRegionBedManifest =
                        serde_json::from_slice(&manifest_bytes).map_err(|error| {
                            region_error(
                                ErrorCode::InvalidInput,
                                format!("could not parse BED manifest: {error}"),
                            )
                        })?;
                    if manifest.schema != gp::GENOMIC_REGION_BED_MANIFEST_SCHEMA {
                        return Err(region_error(
                            ErrorCode::InvalidInput,
                            "unsupported genomic-region BED manifest schema",
                        ));
                    }
                    let expected_columns = BED_COLUMNS
                        .iter()
                        .map(|value| value.to_string())
                        .collect::<Vec<_>>();
                    if manifest.coordinate_convention != "BED6; 0-based, half-open"
                        || manifest.columns != expected_columns
                    {
                        return Err(region_error(
                            ErrorCode::InvalidInput,
                            "BED manifest does not declare the canonical GENtle BED6 column and coordinate contract",
                        ));
                    }
                    if manifest.bed_sha256 != sha256_prefixed_bytes(&bytes) {
                        return Err(region_error(
                            ErrorCode::InvalidInput,
                            "BED content digest does not match its manifest",
                        ));
                    }
                    if manifest.region_set.regions.len() > request.max_rows {
                        return Err(region_error(
                            ErrorCode::InvalidInput,
                            "BED manifest region set exceeds the declared row limit",
                        ));
                    }
                    validate_set_digest(&manifest.region_set)?;
                    if bed_payload(&manifest.region_set).as_bytes() != bytes.as_slice() {
                        return Err(region_error(
                            ErrorCode::InvalidInput,
                            "BED rows do not exactly match the canonical region geometry and labels bound by the manifest",
                        ));
                    }
                    manifest.region_set
                } else {
                    let reference = request.bare_bed_reference.ok_or_else(|| {
                        region_error(
                            ErrorCode::InvalidInput,
                            "bare BED import requires an explicit bound genome reference",
                        )
                    })?;
                    let set_id = request.set_id_override.clone().ok_or_else(|| {
                        region_error(
                            ErrorCode::InvalidInput,
                            "bare BED import requires set_id_override",
                        )
                    })?;
                    Self::parse_bare_bed(&bytes, reference, set_id, request.max_rows)?
                }
            }
        };
        if let Some(set_id) = request.set_id_override {
            valid_identifier(&set_id, "set_id_override")?;
            set.set_id = set_id;
            recompute_set_digest(&mut set)?;
        }
        validate_set_digest(&set)?;
        let mut store = self.genomic_region_store()?;
        if let Some(index) = store.sets.iter().position(|item| item.set_id == set.set_id) {
            if request.collision_policy == gp::GenomicRegionCollisionPolicy::Reject {
                return Err(region_error(
                    ErrorCode::InvalidInput,
                    format!("region set '{}' already exists", set.set_id),
                ));
            }
            store.sets[index] = set.clone();
        } else {
            store.sets.push(set.clone());
        }
        self.persist_genomic_region_store(store)?;
        report_for_region("import", Some(set), None)
    }

    fn export_region_set(
        &self,
        request: gp::GenomicRegionExportRequest,
    ) -> Result<gp::GenomicRegionOperationReport, EngineError> {
        if request.json_path.is_none() && request.bed_path.is_none() {
            return Err(region_error(
                ErrorCode::InvalidInput,
                "export requires json_path and/or bed_path",
            ));
        }
        let store = self.genomic_region_store()?;
        let set = store
            .sets
            .iter()
            .find(|set| set.set_id == request.set_id)
            .cloned()
            .ok_or_else(|| region_error(ErrorCode::NotFound, "genomic region set not found"))?;
        validate_set_digest(&set)?;
        let mut written = vec![];
        if let Some(path) = request.json_path.as_deref() {
            self.write_pretty_json_file(&set, path, "genomic region-set JSON")?;
            written.push(Self::reported_artifact_path(
                path,
                request.include_local_paths,
            ));
        }
        if let Some(path) = request.bed_path.as_deref() {
            let bed = bed_payload(&set);
            self.write_text_file(path, &bed, "genomic region-set BED")?;
            written.push(Self::reported_artifact_path(
                path,
                request.include_local_paths,
            ));
            let manifest_path = request
                .manifest_path
                .clone()
                .unwrap_or_else(|| format!("{path}.manifest.json"));
            let bed_file_name = Path::new(path)
                .file_name()
                .and_then(|value| value.to_str())
                .unwrap_or("regions.bed")
                .to_string();
            let manifest = gp::GenomicRegionBedManifest {
                schema: gp::GENOMIC_REGION_BED_MANIFEST_SCHEMA.to_string(),
                bed_sha256: sha256_prefixed_bytes(bed.as_bytes()),
                bed_file_name,
                coordinate_convention: "BED6; 0-based, half-open".to_string(),
                columns: BED_COLUMNS.iter().map(|value| value.to_string()).collect(),
                region_set: set.clone(),
            };
            self.write_pretty_json_file(
                &manifest,
                &manifest_path,
                "genomic region-set BED manifest",
            )?;
            written.push(Self::reported_artifact_path(
                &manifest_path,
                request.include_local_paths,
            ));
        } else if request.manifest_path.is_some() {
            return Err(region_error(
                ErrorCode::InvalidInput,
                "manifest_path is only valid with bed_path",
            ));
        }
        let mut report = report_for_region("export", Some(set), None)?;
        report.written_artifacts = written;
        Ok(report)
    }

    fn reported_artifact_path(path: &str, include_local_paths: bool) -> String {
        if include_local_paths {
            path.to_string()
        } else {
            Path::new(path)
                .file_name()
                .and_then(|value| value.to_str())
                .unwrap_or("artifact")
                .to_string()
        }
    }

    pub(super) fn apply_genomic_region_operation(
        &mut self,
        op: Operation,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let report = match op {
            Operation::CreateGenomicRegion { request } => self.create_region(request)?,
            Operation::CaptureGenomicRegion { request } => self.capture_region(request)?,
            Operation::ListGenomicRegions { request } => self.list_regions(request)?,
            Operation::InspectGenomicRegion { request } => self.inspect_region(request)?,
            Operation::UpdateGenomicRegionPresentation { request } => {
                self.update_region(request)?
            }
            Operation::DeriveGenomicRegion { request } => self.derive_region(request)?,
            Operation::ImportGenomicRegionSet { request } => self.import_region_set(request)?,
            Operation::ExportGenomicRegionSet { request } => self.export_region_set(request)?,
            _ => {
                return Err(region_error(
                    ErrorCode::Internal,
                    "unexpected operation routed to genomic-region handler",
                ));
            }
        };
        result.messages.push(format!(
            "Genomic-region {} completed{}",
            report.action,
            report
                .region
                .as_ref()
                .map(|region| format!(" for '{}'", region.region_id))
                .unwrap_or_default()
        ));
        result.genomic_region_operation = Some(Box::new(report));
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna_sequence::DNAsequence;
    use tempfile::tempdir;

    fn anchored_engine(anchor_strand: &str) -> GentleEngine {
        let mut engine = GentleEngine::default();
        engine.state_mut().sequences.insert(
            "roi_locus".to_string(),
            DNAsequence::from_sequence(&"ACGT".repeat(25)).expect("synthetic sequence"),
        );
        engine.state_mut().metadata.insert(
            PROVENANCE_METADATA_KEY.to_string(),
            serde_json::json!({
                "genome_extractions": [{
                    "seq_id": "roi_locus",
                    "genome_id": "GRCh38",
                    "chromosome": "chr7",
                    "start_1based": 1001,
                    "end_1based": 1100,
                    "anchor_strand": anchor_strand,
                    "anchor_verified": true,
                    "recorded_at_unix_ms": 123
                }]
            }),
        );
        engine
    }

    fn reference() -> gp::GenomicRegionReference {
        gp::GenomicRegionReference {
            species_scientific_name: Some("Homo sapiens".to_string()),
            taxon_id: Some(9606),
            assembly_name: "GRCh38".to_string(),
            assembly_accession: Some("GCA_000001405.15".to_string()),
            contig_name: "chr7".to_string(),
            contig_accession: Some("NC_000007.14".to_string()),
            contig_aliases: vec!["7".to_string()],
        }
    }

    fn capture_selection(
        engine: &mut GentleEngine,
        set_id: &str,
        start: u64,
        end: u64,
        strand: gp::GenomicRegionStrand,
    ) -> gp::GenomicRegionOperationReport {
        engine
            .apply(Operation::CaptureGenomicRegion {
                request: gp::GenomicRegionCaptureRequest {
                    set_id: set_id.to_string(),
                    set_label: Some("Portable regions".to_string()),
                    source: gp::GenomicRegionCaptureSource::SequenceSelection {
                        seq_id: "roi_locus".to_string(),
                        local_start_0based: start,
                        local_end_0based_exclusive: end,
                        strand,
                        reference_override: Some(reference()),
                    },
                    purpose: gp::GenomicRegionPurpose::CandidateCisRegulatoryRegion,
                    ..Default::default()
                },
            })
            .expect("capture")
            .genomic_region_operation
            .expect("region report")
            .as_ref()
            .clone()
    }

    #[test]
    fn genomic_region_manual_capture_preserves_one_base_and_copy_conventions() {
        let mut engine = anchored_engine("+");
        let report = capture_selection(
            &mut engine,
            "roi_set",
            0,
            1,
            gp::GenomicRegionStrand::Unstranded,
        );
        let region = report.region.expect("region");
        assert_eq!(region.interval.start_0based, 1000);
        assert_eq!(region.interval.end_0based_exclusive, 1001);
        assert_eq!(region.interval.strand, gp::GenomicRegionStrand::Unstranded);
        let expected_bed = format!("chr7\t1000\t1001\t{}\t0\t.", region.region_id);
        assert_eq!(report.bed_row.as_deref(), Some(expected_bed.as_str()));
        let human = report.human_copy.expect("human copy");
        assert!(human.contains("chr7:1001-1001, 1-based inclusive"));
        let round_trip: gp::GenomicRegionOfInterest = serde_json::from_str(
            report
                .canonical_roi_json
                .as_deref()
                .expect("canonical JSON"),
        )
        .expect("canonical ROI JSON");
        assert_eq!(round_trip, region);
        validate_region_digest(&round_trip).expect("valid digest");
    }

    #[test]
    fn genomic_region_reverse_anchor_projects_once_and_keeps_local_strand_explicit() {
        let mut engine = anchored_engine("-");
        let report = capture_selection(
            &mut engine,
            "reverse_set",
            10,
            20,
            gp::GenomicRegionStrand::Plus,
        );
        let region = report.region.expect("region");
        assert_eq!(
            (
                region.interval.start_0based,
                region.interval.end_0based_exclusive
            ),
            (1080, 1090)
        );
        assert_eq!(region.interval.strand, gp::GenomicRegionStrand::Minus);
        let projection = region.local_projection.expect("local projection");
        assert_eq!(
            (
                projection.local_start_0based,
                projection.local_end_0based_exclusive
            ),
            (10, 20)
        );
        assert_eq!(projection.local_strand, gp::GenomicRegionStrand::Plus);
        assert_eq!(projection.anchor_strand, gp::GenomicRegionStrand::Minus);
    }

    #[test]
    fn genomic_region_labels_update_without_changing_identity() {
        let mut engine = anchored_engine("+");
        let captured = capture_selection(
            &mut engine,
            "update_set",
            5,
            15,
            gp::GenomicRegionStrand::Unstranded,
        );
        let original = captured.region.expect("region");
        let updated = engine
            .apply(Operation::UpdateGenomicRegionPresentation {
                request: gp::GenomicRegionUpdateRequest {
                    set_id: "update_set".to_string(),
                    region_id: original.region_id.clone(),
                    label: Some("Readable label".to_string()),
                    display_color_hex: Some("#1a2b3c".to_string()),
                    notes: vec!["reviewed by caller".to_string()],
                    ..Default::default()
                },
            })
            .expect("update")
            .genomic_region_operation
            .expect("report")
            .region
            .expect("updated region");
        assert_eq!(updated.identity_sha256, original.identity_sha256);
        assert_ne!(updated.content_sha256, original.content_sha256);
        assert_eq!(updated.label.as_deref(), Some("Readable label"));
        assert_eq!(updated.display_color_hex.as_deref(), Some("#1A2B3C"));

        let recolored_content_sha256 = updated.content_sha256.clone();
        let preserved = engine
            .apply(Operation::UpdateGenomicRegionPresentation {
                request: gp::GenomicRegionUpdateRequest {
                    set_id: "update_set".to_string(),
                    region_id: updated.region_id.clone(),
                    label: Some("Renamed without a colour field".to_string()),
                    ..Default::default()
                },
            })
            .expect("legacy presentation update")
            .genomic_region_operation
            .expect("report")
            .region
            .expect("updated region");
        assert_eq!(preserved.display_color_hex.as_deref(), Some("#1A2B3C"));
        assert_ne!(preserved.content_sha256, recolored_content_sha256);

        let invalid = engine.apply(Operation::UpdateGenomicRegionPresentation {
            request: gp::GenomicRegionUpdateRequest {
                set_id: "update_set".to_string(),
                region_id: updated.region_id,
                display_color_hex: Some("orange".to_string()),
                ..Default::default()
            },
        });
        assert!(invalid.is_err());
    }

    #[test]
    fn genomic_region_derivations_are_explicit_and_disjoint_union_is_rejected() {
        let mut engine = anchored_engine("+");
        let first = capture_selection(
            &mut engine,
            "derive_set",
            10,
            30,
            gp::GenomicRegionStrand::Plus,
        )
        .region
        .expect("first");
        let second = capture_selection(
            &mut engine,
            "derive_set",
            20,
            40,
            gp::GenomicRegionStrand::Plus,
        )
        .region
        .expect("second");
        let intersection = engine
            .apply(Operation::DeriveGenomicRegion {
                request: gp::GenomicRegionDeriveRequest {
                    set_id: "derive_set".to_string(),
                    parent_region_ids: vec![first.region_id.clone(), second.region_id.clone()],
                    method: gp::GenomicRegionDerivationMethod::Intersection,
                    region_id: Some("intersection".to_string()),
                    ..Default::default()
                },
            })
            .expect("intersection")
            .genomic_region_operation
            .expect("report")
            .region
            .expect("region");
        assert_eq!(
            (
                intersection.interval.start_0based,
                intersection.interval.end_0based_exclusive
            ),
            (1020, 1030)
        );
        assert_eq!(
            intersection.derivation.expect("derivation").method,
            gp::GenomicRegionDerivationMethod::Intersection
        );
        let union = engine
            .apply(Operation::DeriveGenomicRegion {
                request: gp::GenomicRegionDeriveRequest {
                    set_id: "derive_set".to_string(),
                    parent_region_ids: vec![first.region_id.clone(), second.region_id],
                    method: gp::GenomicRegionDerivationMethod::Union,
                    region_id: Some("union".to_string()),
                    ..Default::default()
                },
            })
            .expect("contiguous union")
            .genomic_region_operation
            .expect("report")
            .region
            .expect("region");
        assert_eq!(
            (
                union.interval.start_0based,
                union.interval.end_0based_exclusive
            ),
            (1010, 1040)
        );

        let third = capture_selection(
            &mut engine,
            "derive_set",
            60,
            70,
            gp::GenomicRegionStrand::Plus,
        )
        .region
        .expect("third");
        let error = engine
            .apply(Operation::DeriveGenomicRegion {
                request: gp::GenomicRegionDeriveRequest {
                    set_id: "derive_set".to_string(),
                    parent_region_ids: vec![first.region_id.clone(), third.region_id.clone()],
                    method: gp::GenomicRegionDerivationMethod::Union,
                    ..Default::default()
                },
            })
            .expect_err("disjoint union must not silently become a hull");
        assert_eq!(error.code, ErrorCode::InvalidInput);
        assert!(error.message.contains("explicit hull"));
        let hull = engine
            .apply(Operation::DeriveGenomicRegion {
                request: gp::GenomicRegionDeriveRequest {
                    set_id: "derive_set".to_string(),
                    parent_region_ids: vec![first.region_id, third.region_id],
                    method: gp::GenomicRegionDerivationMethod::Hull,
                    region_id: Some("explicit_hull".to_string()),
                    ..Default::default()
                },
            })
            .expect("explicit hull")
            .genomic_region_operation
            .expect("report")
            .region
            .expect("region");
        assert_eq!(
            (
                hull.interval.start_0based,
                hull.interval.end_0based_exclusive
            ),
            (1010, 1070)
        );
        assert!(
            hull.derivation
                .expect("derivation")
                .rule
                .contains("bases not present")
        );
    }

    #[test]
    fn genomic_region_json_and_bed_manifest_round_trip_and_bare_bed_requires_reference() {
        let temp = tempdir().expect("tempdir");
        let json_path = temp.path().join("regions.json");
        let bed_path = temp.path().join("regions.bed");
        let manifest_path = temp.path().join("regions.bed.manifest.json");
        let mut source = anchored_engine("+");
        let captured = capture_selection(
            &mut source,
            "portable_set",
            1,
            9,
            gp::GenomicRegionStrand::Plus,
        );
        let region_id = captured.region.expect("region").region_id;
        source
            .apply(Operation::UpdateGenomicRegionPresentation {
                request: gp::GenomicRegionUpdateRequest {
                    set_id: "portable_set".to_string(),
                    region_id,
                    label: Some("unsafe label / tabs\tbecome safe".to_string()),
                    ..Default::default()
                },
            })
            .expect("label update");
        let export = source
            .apply(Operation::ExportGenomicRegionSet {
                request: gp::GenomicRegionExportRequest {
                    set_id: "portable_set".to_string(),
                    json_path: Some(json_path.display().to_string()),
                    bed_path: Some(bed_path.display().to_string()),
                    manifest_path: Some(manifest_path.display().to_string()),
                    include_local_paths: false,
                },
            })
            .expect("export")
            .genomic_region_operation
            .expect("report");
        assert_eq!(
            export.written_artifacts,
            vec![
                "regions.json".to_string(),
                "regions.bed".to_string(),
                "regions.bed.manifest.json".to_string()
            ]
        );
        let bed = std::fs::read_to_string(&bed_path).expect("BED");
        assert!(bed.contains("unsafe_label___tabs_become_safe"));

        let mut json_target = GentleEngine::default();
        json_target
            .apply(Operation::ImportGenomicRegionSet {
                request: gp::GenomicRegionImportRequest {
                    path: json_path.display().to_string(),
                    format: gp::GenomicRegionImportFormat::Json,
                    ..Default::default()
                },
            })
            .expect("JSON import");
        assert_eq!(
            json_target
                .genomic_region_store_snapshot()
                .expect("store")
                .sets[0],
            source
                .genomic_region_store_snapshot()
                .expect("source store")
                .sets[0]
        );

        let mut bed_target = GentleEngine::default();
        bed_target
            .apply(Operation::ImportGenomicRegionSet {
                request: gp::GenomicRegionImportRequest {
                    path: bed_path.display().to_string(),
                    format: gp::GenomicRegionImportFormat::Bed,
                    manifest_path: Some(manifest_path.display().to_string()),
                    ..Default::default()
                },
            })
            .expect("manifest-bound BED import");
        let altered_bed_path = temp.path().join("regions.altered.bed");
        std::fs::write(&altered_bed_path, bed.replacen("\t1001\t", "\t1002\t", 1))
            .expect("write altered BED");
        let manifest_error = GentleEngine::default()
            .apply(Operation::ImportGenomicRegionSet {
                request: gp::GenomicRegionImportRequest {
                    path: altered_bed_path.display().to_string(),
                    format: gp::GenomicRegionImportFormat::Bed,
                    manifest_path: Some(manifest_path.display().to_string()),
                    ..Default::default()
                },
            })
            .expect_err("altered BED must not pass the original manifest");
        assert!(manifest_error.message.contains("digest"));
        let bare_error = GentleEngine::default()
            .apply(Operation::ImportGenomicRegionSet {
                request: gp::GenomicRegionImportRequest {
                    path: bed_path.display().to_string(),
                    format: gp::GenomicRegionImportFormat::Bed,
                    set_id_override: Some("bare".to_string()),
                    ..Default::default()
                },
            })
            .expect_err("bare BED needs an explicit reference");
        assert!(
            bare_error
                .message
                .contains("explicit bound genome reference")
        );
    }

    #[test]
    fn genomic_region_source_capture_retains_conservative_cutrun_and_ensembl_evidence() {
        let mut engine = anchored_engine("+");
        let cutrun_report = gp::CutRunRegulatorySupportReport {
            schema: "gentle.cutrun_regulatory_support.v1".to_string(),
            seq_id: "roi_locus".to_string(),
            evidence_sources: vec![
                gp::CutRunRegulatoryEvidenceSourceRef {
                    source_id: "synthetic_cutrun_dataset".to_string(),
                    dataset_id: Some("dataset_1".to_string()),
                    target_factor: Some("TF-DEMO".to_string()),
                    description: Some("synthetic occupancy lane".to_string()),
                    ..Default::default()
                },
                gp::CutRunRegulatoryEvidenceSourceRef {
                    source_kind: gp::CutRunRegulatoryEvidenceSourceKind::ReadReport,
                    source_id: "synthetic_cutrun_reads".to_string(),
                    report_id: Some("cutrun_report_1".to_string()),
                    target_factor: Some("TF-DEMO".to_string()),
                    ..Default::default()
                },
            ],
            support_windows: vec![gp::CutRunSupportWindowRecord {
                window_id: "support_1".to_string(),
                local_start_0based: 10,
                local_end_0based_exclusive: 20,
                genomic_start_1based: 1011,
                genomic_end_1based: 1020,
                support_strength: gp::CutRunSupportStrength::Strong,
                max_signal_value: Some(4.5),
                ..Default::default()
            }],
            ..Default::default()
        };
        let cutrun = engine
            .apply(Operation::CaptureGenomicRegion {
                request: gp::GenomicRegionCaptureRequest {
                    set_id: "evidence_set".to_string(),
                    source: gp::GenomicRegionCaptureSource::CutrunSupportWindow {
                        report: Box::new(cutrun_report),
                        window_id: "support_1".to_string(),
                        reference: reference(),
                    },
                    purpose: gp::GenomicRegionPurpose::OccupancyRegion,
                    ..Default::default()
                },
            })
            .expect("CUT&RUN capture")
            .genomic_region_operation
            .expect("report")
            .region
            .expect("region");
        let cutrun_evidence = &cutrun.evidence[0];
        assert_eq!(
            cutrun.selection_method,
            gp::GenomicRegionSelectionMethod::CutrunSupportWindow
        );
        assert_eq!(cutrun.evidence.len(), 2);
        assert_eq!(cutrun.evidence[0].source_id, "synthetic_cutrun_dataset");
        assert_eq!(cutrun.evidence[1].source_id, "synthetic_cutrun_reads");
        assert_eq!(
            cutrun.evidence[1].report_id.as_deref(),
            Some("cutrun_report_1")
        );
        assert!(
            cutrun
                .evidence
                .iter()
                .all(|evidence| evidence.assay.as_deref() == Some("CUT&RUN"))
        );
        assert!(
            cutrun_evidence
                .evidence_statement
                .contains("occupancy evidence")
        );
        assert!(
            cutrun_evidence
                .non_claims
                .iter()
                .any(|claim| claim.contains("not biochemical affinity"))
        );
        assert!(
            cutrun_evidence
                .non_claims
                .iter()
                .any(|claim| claim.contains("causal"))
        );

        let source = gp::EnsemblRegulationSourceDescriptor {
            source_id: "synthetic_ensembl_regulation".to_string(),
            provider: "Ensembl Regulation".to_string(),
            annotation_release: "synthetic-r1".to_string(),
            species_scientific_name: "Homo sapiens".to_string(),
            taxon_id: 9606,
            assembly_name: "GRCh38".to_string(),
            assembly_accession: "GCA_000001405.15".to_string(),
            feature_page_url_template:
                "https://regulation.ensembl.org/homo_sapiens/RegulatoryFeature/Summary?rf={FEATURE_ID}"
                    .to_string(),
            ..Default::default()
        };
        let ensembl_report = gp::EnsemblRegulationOverlapReport {
            schema: gp::ENSEMBL_REGULATION_OVERLAP_SCHEMA.to_string(),
            report_id: "ensembl_overlap_1".to_string(),
            seq_id: "roi_locus".to_string(),
            source,
            content_identity_verified: true,
            rows: vec![gp::EnsemblRegulationOverlapRow {
                interval: gp::EnsemblRegulationInterval {
                    chromosome: "chr7".to_string(),
                    start_0based: 1030,
                    end_0based_exclusive: 1040,
                    feature_id: "ENSR_DEMO_1".to_string(),
                    feature_type: "promoter".to_string(),
                    associated_gene_names: vec!["DEMO".to_string()],
                    ..Default::default()
                },
                local_start_0based: 30,
                local_end_0based_exclusive: 40,
                genomic_start_0based: 1030,
                genomic_end_0based_exclusive: 1040,
                overlap_bp: 10,
                ..Default::default()
            }],
            evidence_statement: "Provider annotation overlaps the requested locus.".to_string(),
            non_claims: vec![
                "Associated genes are provider annotations, not causal assignments.".to_string(),
            ],
            ..Default::default()
        };
        let ensembl = engine
            .apply(Operation::CaptureGenomicRegion {
                request: gp::GenomicRegionCaptureRequest {
                    set_id: "evidence_set".to_string(),
                    source: gp::GenomicRegionCaptureSource::EnsemblRegulatoryFeature {
                        report: Box::new(ensembl_report),
                        feature_id: "ENSR_DEMO_1".to_string(),
                        interval_kind: gp::GenomicRegionEnsemblIntervalKind::Core,
                    },
                    ..Default::default()
                },
            })
            .expect("Ensembl capture")
            .genomic_region_operation
            .expect("report")
            .region
            .expect("region");
        let evidence = &ensembl.evidence[0];
        assert_eq!(
            evidence.feature_or_window_id.as_deref(),
            Some("ENSR_DEMO_1")
        );
        assert_eq!(evidence.associated_gene_names, vec!["DEMO"]);
        assert!(
            evidence
                .provider_url
                .as_deref()
                .is_some_and(|url| url.contains("ENSR_DEMO_1"))
        );
        assert!(evidence.non_claims[0].contains("not causal"));
    }

    #[test]
    fn genomic_region_import_validation_rejects_signed_projection_and_reference_mismatches() {
        let mut engine = anchored_engine("+");
        let first = capture_selection(
            &mut engine,
            "validated_set",
            10,
            20,
            gp::GenomicRegionStrand::Plus,
        )
        .region
        .expect("first region");

        let mut inconsistent_projection = first.clone();
        inconsistent_projection
            .local_projection
            .as_mut()
            .expect("projection")
            .local_start_0based += 1;
        recompute_region_digests(&mut inconsistent_projection).expect("signed fixture");
        let error = validate_region_digest(&inconsistent_projection)
            .expect_err("a digest must not legitimize inconsistent projection geometry");
        assert!(error.message.contains("does not reproduce"));

        let mut other_assembly = first;
        other_assembly.region_id = "other_assembly".to_string();
        other_assembly.interval.reference.assembly_name = "T2T-CHM13v2.0".to_string();
        other_assembly.local_projection = None;
        recompute_region_digests(&mut other_assembly).expect("other assembly fixture");
        let mut mixed = engine
            .genomic_region_store_snapshot()
            .expect("store")
            .sets
            .remove(0);
        mixed.regions.push(other_assembly);
        recompute_set_digest(&mut mixed).expect("mixed set digest");
        let error = validate_set_digest(&mixed)
            .expect_err("a set must not mix assemblies even when all region digests are valid");
        assert!(error.message.contains("mixes incompatible"));

        let mut duplicate = engine
            .genomic_region_store_snapshot()
            .expect("store")
            .sets
            .remove(0);
        duplicate.regions.push(duplicate.regions[0].clone());
        recompute_set_digest(&mut duplicate).expect("duplicate set digest");
        let error = validate_set_digest(&duplicate)
            .expect_err("duplicate region identifiers must be rejected");
        assert!(error.message.contains("duplicate region_id"));
    }

    #[test]
    fn genomic_region_stale_evidence_remains_inspectable() {
        let mut engine = GentleEngine::default();
        let report = engine
            .apply(Operation::CaptureGenomicRegion {
                request: gp::GenomicRegionCaptureRequest {
                    set_id: "stale_set".to_string(),
                    source: gp::GenomicRegionCaptureSource::ProviderAnnotation {
                        interval: gp::GenomicRegionInterval {
                            reference: reference(),
                            start_0based: 0,
                            end_0based_exclusive: 1,
                            ..Default::default()
                        },
                        evidence: gp::GenomicRegionEvidenceReference {
                            evidence_id: "provider:stale".to_string(),
                            source_kind: "provider_annotation".to_string(),
                            source_id: "provider_snapshot".to_string(),
                            availability: gp::GenomicRegionEvidenceAvailability::Stale,
                            evidence_statement:
                                "The original provider snapshot is unavailable for revalidation."
                                    .to_string(),
                            ..Default::default()
                        },
                    },
                    ..Default::default()
                },
            })
            .expect("capture stale evidence")
            .genomic_region_operation
            .expect("report");
        let region_id = report.region.as_ref().expect("region").region_id.clone();
        assert_eq!(
            report
                .region
                .as_ref()
                .expect("region")
                .interval
                .start_0based,
            0
        );
        assert!(
            report
                .human_copy
                .as_deref()
                .is_some_and(|copy| copy.contains("chr7:1-1, 1-based inclusive"))
        );
        let inspected = engine
            .apply(Operation::InspectGenomicRegion {
                request: gp::GenomicRegionInspectRequest {
                    set_id: "stale_set".to_string(),
                    region_id,
                },
            })
            .expect("inspect")
            .genomic_region_operation
            .expect("report")
            .region
            .expect("region");
        assert_eq!(
            inspected.evidence[0].availability,
            gp::GenomicRegionEvidenceAvailability::Stale
        );
    }
}
