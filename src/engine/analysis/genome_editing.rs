//! Strict, read-only SNV intent validation over existing reference primitives.
//! Coordinates describe a proposed allele, never an editor, genotype or outcome.

use super::*;
use gentle_protocol::genome_editing::*;
use gentle_protocol::{BiologicalContext, GenomicRegionInterval, GenomicRegionStrand};

fn refusal(reason: &str, detail: impl Into<String>) -> EngineError {
    EngineError {
        code: ErrorCode::InvalidInput,
        message: detail.into(),
        cause_chain: vec![format!("genome_edit_intent.{reason}")],
    }
}

fn text(value: &str, field: &str) -> Result<(), EngineError> {
    if value.is_empty()
        || value.len() > MAX_GENOME_EDIT_ID_BYTES
        || value.trim() != value
        || value.chars().any(char::is_control)
    {
        return Err(refusal(
            "invalid_identity",
            format!("{field} requires 1..256 trimmed, non-control bytes"),
        ));
    }
    Ok(())
}

fn versioned_accession(value: Option<&str>, field: &str) -> Result<(), EngineError> {
    let value = value
        .ok_or_else(|| refusal("missing_reference_identity", format!("{field} is required")))?;
    text(value, field)?;
    if !value.rsplit_once('.').is_some_and(|(id, version)| {
        !id.is_empty()
            && id.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_')
            && !version.is_empty()
            && version.bytes().all(|b| b.is_ascii_digit())
            && version.parse::<u64>().is_ok_and(|v| v > 0)
    }) {
        return Err(refusal(
            "unversioned_reference",
            format!("{field} requires an explicit accession.version"),
        ));
    }
    Ok(())
}

fn dna(value: &str, field: &str) -> Result<String, EngineError> {
    if !value
        .bytes()
        .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    {
        return Err(refusal(
            "ambiguous_dna",
            format!("{field} must contain only A/C/G/T; no cleaning, ambiguity or U conversion"),
        ));
    }
    Ok(value.to_ascii_uppercase())
}

fn context_status(
    reference: &GenomeEditReference,
    context: Option<&BiologicalContext>,
) -> Result<GenomeEditAnnotationStatus, EngineError> {
    let Some(context) = context else {
        return Ok(GenomeEditAnnotationStatus::AnnotationNotSupplied);
    };
    text(&context.context_id, "context.context_id")?;
    for (field, value) in [
        ("organism", &context.organism),
        ("taxon_id", &context.taxon_id),
        ("genome_id", &context.genome_id),
        ("assembly_accession", &context.assembly_accession),
        ("annotation_source", &context.annotation_source),
        ("annotation_release", &context.annotation_release),
        ("symbol_namespace", &context.symbol_namespace),
    ] {
        if let Some(value) = value {
            text(value, field)?;
        }
    }
    match reference {
        GenomeEditReference::Genomic { reference } => {
            if context
                .assembly_accession
                .as_ref()
                .is_some_and(|v| Some(v) != reference.assembly_accession.as_ref())
            {
                return Err(refusal(
                    "context_conflict",
                    "Context assembly differs from the exact reference",
                ));
            }
            if let Some(taxon) = &context.taxon_id {
                let parsed = taxon.parse::<u32>().map_err(|_| {
                    refusal(
                        "context_conflict",
                        "Context taxon must be a numeric u32 identity",
                    )
                })?;
                if Some(parsed) != reference.taxon_id {
                    return Err(refusal(
                        "context_conflict",
                        "Context taxon differs from the exact reference",
                    ));
                }
            }
        }
        GenomeEditReference::Synthetic { .. } => {
            if context.organism.is_some()
                || context.taxon_id.is_some()
                || context.assembly_accession.is_some()
                || context.genome_id.is_some()
            {
                return Err(refusal(
                    "synthetic_identity_conflict",
                    "Synthetic reference cannot carry biological organism/taxon/assembly/genome identities",
                ));
            }
        }
    }
    match (&context.annotation_source, &context.annotation_release) {
        (None, None) => Ok(GenomeEditAnnotationStatus::AnnotationNotSupplied),
        (Some(_), Some(_)) => Ok(GenomeEditAnnotationStatus::SuppliedNotEvaluated),
        _ => Err(refusal(
            "incomplete_annotation_binding",
            "Bound annotations require both source and release; omit both for sequence-only validation",
        )),
    }
}

impl GentleEngine {
    /// Validate one bounded SNV without files, providers, persistence or history.
    pub fn evaluate_genome_edit_intent(
        &self,
        mut request: GenomeEditIntentRequest,
    ) -> Result<GenomeEditIntentReport, EngineError> {
        if request.schema != GENOME_EDIT_INTENT_REQUEST_SCHEMA {
            return Err(refusal(
                "schema",
                "Expected gentle.genome_edit_intent_request.v1",
            ));
        }
        text(&request.design_haplotype, "design_haplotype")?;
        if ["unknown", "unresolved", "?"]
            .iter()
            .any(|s| request.design_haplotype.eq_ignore_ascii_case(s))
        {
            return Err(refusal(
                "unresolved_haplotype",
                "Select an explicit design haplotype, not an unresolved sample genotype",
            ));
        }
        match &request.reference {
            GenomeEditReference::Genomic { reference } => {
                text(&reference.assembly_name, "assembly_name")?;
                text(&reference.contig_name, "contig_name")?;
                versioned_accession(
                    reference.assembly_accession.as_deref(),
                    "assembly_accession",
                )?;
                versioned_accession(reference.contig_accession.as_deref(), "contig_accession")?;
                if !reference.taxon_id.is_some_and(|id| id > 0) {
                    return Err(refusal(
                        "missing_reference_identity",
                        "Genomic reference requires a positive taxon_id",
                    ));
                }
                if let Some(species) = &reference.species_scientific_name {
                    text(species, "species_scientific_name")?;
                }
                if reference.contig_aliases.len() > 32 {
                    return Err(refusal(
                        "identity_budget",
                        "At most 32 explicit contig aliases are allowed",
                    ));
                }
                for alias in &reference.contig_aliases {
                    text(alias, "contig alias")?;
                }
            }
            GenomeEditReference::Synthetic {
                namespace,
                reference_id,
            } => {
                text(namespace, "namespace")?;
                text(reference_id, "reference_id")?;
                if !namespace
                    .strip_prefix("synthetic:")
                    .is_some_and(|s| !s.is_empty())
                {
                    return Err(refusal(
                        "synthetic_namespace",
                        "Synthetic reference requires namespace=synthetic:<explicit-id>",
                    ));
                }
            }
        }
        let annotation_status = context_status(&request.reference, request.context.as_ref())?;
        let hex = request
            .expected_sequence_sha256
            .strip_prefix("sha256:")
            .ok_or_else(|| {
                refusal(
                    "digest_format",
                    "Expected sha256: followed by 64 hexadecimal digits",
                )
            })?;
        if hex.len() != 64 || !hex.bytes().all(|b| b.is_ascii_hexdigit()) {
            return Err(refusal(
                "digest_format",
                "Expected sha256: followed by 64 hexadecimal digits",
            ));
        }
        let anchor = &request.anchor;
        if anchor.start_1based == 0
            || anchor.start_1based > anchor.end_1based
            || anchor.strand == GenomicRegionStrand::Unstranded
        {
            return Err(refusal(
                "anchor_mismatch",
                "Anchor must be positive 1-based inclusive with explicit plus/minus strand",
            ));
        }
        let origin = anchor
            .start_1based
            .checked_sub(1)
            .ok_or_else(|| refusal("anchor_mismatch", "Anchor origin underflows"))?;
        let anchor_len = anchor
            .end_1based
            .checked_sub(anchor.start_1based)
            .and_then(|n| n.checked_add(1))
            .ok_or_else(|| refusal("anchor_mismatch", "Anchor length overflows"))?;
        if anchor_len > MAX_GENOME_EDIT_WINDOW_BP as u64 {
            return Err(refusal(
                "window_budget",
                "SNV validation admits windows of at most 100000 bases",
            ));
        }
        let source_sequence = match &request.target {
            GenomeEditTarget::InlineSequence { sequence } => {
                if sequence.is_empty() || sequence.len() > MAX_GENOME_EDIT_WINDOW_BP {
                    return Err(refusal(
                        "window_budget",
                        "Inline DNA requires 1..100000 bases",
                    ));
                }
                dna(sequence, "inline sequence")?
            }
            GenomeEditTarget::LoadedSequence { seq_id } => {
                text(seq_id, "seq_id")?;
                let sequence = self
                    .state
                    .sequences
                    .get(seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                        cause_chain: vec!["genome_edit_intent.sequence_unavailable".into()],
                    })?;
                if sequence.is_circular()
                    || sequence
                        .molecule_type()
                        .is_some_and(|s| !s.eq_ignore_ascii_case("DNA"))
                {
                    return Err(EngineError::new(
                        ErrorCode::Unsupported,
                        "SNV intent validation requires linear DNA",
                    ));
                }
                if sequence.len() == 0 || sequence.len() > MAX_GENOME_EDIT_WINDOW_BP {
                    return Err(refusal(
                        "window_budget",
                        "Loaded DNA requires 1..100000 bases",
                    ));
                }
                dna(&sequence.get_forward_string(), "loaded sequence")?
            }
        };
        if u64::try_from(source_sequence.len()).ok() != Some(anchor_len) {
            return Err(refusal(
                "anchor_mismatch",
                "Inclusive anchor length does not equal the exact source sequence length",
            ));
        }
        let source_sequence_sha256 = sha256_prefixed_str(&source_sequence);
        if !source_sequence_sha256.eq_ignore_ascii_case(&request.expected_sequence_sha256) {
            return Err(refusal(
                "sequence_digest_mismatch",
                "Expected source sequence digest does not match current uppercase DNA",
            ));
        }
        if request.intent.deleted.len() > MAX_GENOME_EDIT_WINDOW_BP
            || request.intent.inserted.len() > MAX_GENOME_EDIT_WINDOW_BP
        {
            return Err(refusal(
                "window_budget",
                "Alleles exceed the 100000-base admission bound",
            ));
        }
        request.intent.deleted = dna(&request.intent.deleted, "deleted allele")?;
        request.intent.inserted = dna(&request.intent.inserted, "inserted allele")?;
        let edit = &request.intent;
        let span = edit
            .end_0based_exclusive
            .checked_sub(edit.start_0based)
            .ok_or_else(|| refusal("interval", "Edit interval is reversed"))?;
        if span != edit.deleted.len() as u64 {
            return Err(refusal(
                "allele_span",
                "Deleted allele length must equal the genomic half-open span",
            ));
        }
        if edit.deleted == edit.inserted {
            return Err(refusal(
                "unchanged_intent",
                "Intent must change the reference allele",
            ));
        }
        if edit.start_0based < origin || edit.end_0based_exclusive > anchor.end_1based {
            return Err(refusal(
                "outside_sequence",
                "Edit lies outside the exact source window",
            ));
        }
        if span != 1 || edit.inserted.len() != 1 {
            return Err(EngineError::new(
                ErrorCode::Unsupported,
                "S1 admits one SNV only; indels and multi-base edits require a later contract",
            ));
        }
        let genomic_start = usize::try_from(edit.start_0based - origin)
            .map_err(|_| refusal("index_overflow", "Local start does not fit usize"))?;
        let genomic_end = usize::try_from(edit.end_0based_exclusive - origin)
            .map_err(|_| refusal("index_overflow", "Local end does not fit usize"))?;
        let reversed = anchor.strand == GenomicRegionStrand::Minus;
        let genomic_forward_sequence = if reversed {
            Self::reverse_complement(&source_sequence)
        } else {
            source_sequence.clone()
        };
        if genomic_forward_sequence.get(genomic_start..genomic_end) != Some(edit.deleted.as_str()) {
            return Err(refusal(
                "reference_base_mismatch",
                "Deleted genomic-forward base does not match the bound reference",
            ));
        }
        let (local_start, local_end) = if reversed {
            (
                source_sequence.len() - genomic_end,
                source_sequence.len() - genomic_start,
            )
        } else {
            (genomic_start, genomic_end)
        };
        let mut local_projection = None;
        let reference_verification = match (&request.target, &request.reference) {
            (
                GenomeEditTarget::LoadedSequence { seq_id },
                GenomeEditReference::Genomic { reference },
            ) => {
                // Typed provenance refuses overflowing legacy usize anchors before
                // the older summary/projection helpers can cast JSON coordinates.
                let provenance = self
                    .latest_genome_extraction_provenance_for_seq(seq_id)
                    .ok_or_else(|| {
                        refusal(
                            "anchor_mismatch",
                            "Loaded genomic DNA requires complete, typed retained provenance",
                        )
                    })?;
                if Some(&provenance.genome_id) != reference.assembly_accession.as_ref()
                    || provenance.chromosome.as_ref() != reference.contig_accession.as_ref()
                    || provenance.start_1based.and_then(|v| u64::try_from(v).ok())
                        != Some(anchor.start_1based)
                    || provenance.end_1based.and_then(|v| u64::try_from(v).ok())
                        != Some(anchor.end_1based)
                    || provenance.anchor_strand != Some(if reversed { '-' } else { '+' })
                    || provenance.anchor_verified != Some(true)
                {
                    return Err(refusal(
                        "anchor_mismatch",
                        "Latest typed provenance must match the exact versioned assembly, contig, coordinates, strand and verified flag; no older-record or alias fallback",
                    ));
                }
                let retained = self.sequence_genome_anchor_summary(seq_id).map_err(|_| {
                    refusal(
                        "anchor_mismatch",
                        "Loaded genomic DNA requires retained exact anchor provenance",
                    )
                })?;
                if retained.anchor_verified != Some(true)
                    || Some(&retained.genome_id) != reference.assembly_accession.as_ref()
                    || Some(&retained.chromosome) != reference.contig_accession.as_ref()
                    || u64::try_from(retained.start_1based).ok() != Some(anchor.start_1based)
                    || u64::try_from(retained.end_1based).ok() != Some(anchor.end_1based)
                    || retained.strand != Some(if reversed { '-' } else { '+' })
                {
                    return Err(refusal(
                        "anchor_mismatch",
                        "Requested anchor differs from verified retained provenance",
                    ));
                }
                if request
                    .context
                    .as_ref()
                    .and_then(|c| c.genome_id.as_ref())
                    .is_some_and(|id| id != &retained.genome_id)
                {
                    return Err(refusal(
                        "context_conflict",
                        "Context genome_id differs from retained provenance",
                    ));
                }
                let interval = GenomicRegionInterval {
                    reference: reference.clone(),
                    start_0based: edit.start_0based,
                    end_0based_exclusive: edit.end_0based_exclusive,
                    strand: GenomicRegionStrand::Plus,
                    ..Default::default()
                };
                let projection = self
                    .local_projection_for_interval(seq_id, &interval)
                    .map_err(|error| refusal("anchor_mismatch", error.message))?;
                // The existing projection binds stored bytes, while this
                // request binds uppercase DNA without modifying those bytes.
                let projection_source_sha256 = self
                    .state
                    .sequences
                    .get(seq_id)
                    .map(|dna| sha256_prefixed_bytes(dna.forward_bytes()));
                if projection.status != gentle_protocol::GenomicRegionLocalProjectionStatus::Current
                    || projection_source_sha256.as_deref()
                        != Some(projection.sequence_sha256.as_str())
                    || projection.local_start_0based != local_start as u64
                    || projection.local_end_0based_exclusive != local_end as u64
                {
                    return Err(refusal(
                        "anchor_mismatch",
                        "Engine projection failed fresh digest/geometry checks",
                    ));
                }
                local_projection = Some(projection);
                GenomeEditReferenceVerification::RetainedAnchorMatched
            }
            (
                GenomeEditTarget::LoadedSequence { seq_id },
                GenomeEditReference::Synthetic { .. },
            ) => {
                if self.sequence_genome_anchor_summary(seq_id).is_ok() {
                    return Err(refusal(
                        "synthetic_identity_conflict",
                        "Genome-anchored DNA cannot be silently relabelled synthetic",
                    ));
                }
                GenomeEditReferenceVerification::SyntheticDeclared
            }
            (_, GenomeEditReference::Synthetic { .. }) => {
                GenomeEditReferenceVerification::SyntheticDeclared
            }
            (_, GenomeEditReference::Genomic { .. }) => {
                GenomeEditReferenceVerification::InlineCallerDeclared
            }
        };
        let desired_genomic_forward_sequence = format!(
            "{}{}{}",
            &genomic_forward_sequence[..genomic_start],
            edit.inserted,
            &genomic_forward_sequence[genomic_end..]
        );
        let desired_local_sequence = if reversed {
            Self::reverse_complement(&desired_genomic_forward_sequence)
        } else {
            desired_genomic_forward_sequence.clone()
        };
        Ok(GenomeEditIntentReport {
            schema: GENOME_EDIT_INTENT_REPORT_SCHEMA.into(), request,
            source_sequence_sha256,
            genomic_forward_sha256: sha256_prefixed_str(&genomic_forward_sequence),
            desired_genomic_forward_sha256: sha256_prefixed_str(&desired_genomic_forward_sequence),
            desired_local_sha256: sha256_prefixed_str(&desired_local_sequence),
            source_sequence, genomic_forward_sequence, desired_genomic_forward_sequence, desired_local_sequence,
            local_start_0based: local_start as u64, local_end_0based_exclusive: local_end as u64,
            local_projection, reference_verification, annotation_status,
            nonclaims: vec![
                "Sequence/coordinate consistency only; no editor, guide, activity window, efficiency or off-target assessment.".into(),
                "Reference metadata and retained anchors are not independently authenticated against a provider; no network or annotation query was performed.".into(),
                "The selected haplotype is a design basis, not a measured cell-line genotype. Annotation labels are not evaluated consequences.".into(),
                "The desired allele is returned as data, not created in the project or observed experimentally; no laboratory or functional approval.".into(),
            ],
        })
    }
}

#[cfg(test)]
#[path = "genome_editing/tests.rs"]
mod tests;
