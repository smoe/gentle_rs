//! Hash-bound locus/TATA evidence projected onto already scored TSS windows.
//!
//! No motif scoring, peak calling, resource retrieval or project mutation occurs.

use super::{GentleEngine, invalid};
use crate::{digest_utils::sha256_hex_bytes, engine::EngineError, locus_report::LocusDocument};
use gentle_protocol::{isoform_evidence::*, tata_boxes::*, tss_profiles::*};
use std::{collections::BTreeSet, io::Read, path::Path};

const MAX_INPUT_BYTES: u64 = 128 * 1024 * 1024;

#[cfg(test)]
#[path = "tss_profiles_context_tests.rs"]
mod tests;

fn read(path: &Path, limit: u64) -> Result<Vec<u8>, EngineError> {
    let file = crate::tss_fasta_bundle::open_regular_input(path)?;
    let mut bytes = Vec::new();
    file.take(limit + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| invalid(format!("TSS context {}: {e}", path.display())))?;
    if bytes.len() as u64 > limit {
        return Err(invalid("TSS context input exceeds its bounded file limit"));
    }
    Ok(bytes)
}

fn digest(raw: &str) -> Result<&str, EngineError> {
    let value = raw.strip_prefix("sha256:").unwrap_or(raw);
    gentle_engine::tss_profiles::validate_sha256(value)?;
    Ok(value)
}

fn bound_file(base: &Path, input: &TssContextFile) -> Result<Vec<u8>, EngineError> {
    if input.path.trim().is_empty() {
        return Err(invalid("TSS context requires a nonempty input path"));
    }
    let bytes = read(&base.join(&input.path), MAX_INPUT_BYTES)?;
    if sha256_hex_bytes(&bytes) != digest(&input.sha256)? {
        return Err(invalid(format!(
            "TSS context input SHA-256 mismatch: {}",
            input.path
        )));
    }
    Ok(bytes)
}

fn chromosome(value: &str) -> &str {
    value.strip_prefix("chr").unwrap_or(value)
}

fn local_genomic(
    anchor: &GeneLocusGenomeAnchorBinding,
    start: usize,
    end: usize,
) -> Result<(u64, u64), EngineError> {
    let length = anchor
        .end_1based
        .checked_sub(anchor.start_1based)
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| invalid("Invalid locus anchor"))?;
    if start == 0 || end < start || end > length {
        return Err(invalid(
            "TSS context source interval is outside the bound locus sequence",
        ));
    }
    let pair = match anchor.strand {
        Some('+') => (
            anchor.start_1based + (start - 1),
            anchor.start_1based + (end - 1),
        ),
        Some('-') => (anchor.end_1based - end + 1, anchor.end_1based - start + 1),
        _ => {
            return Err(invalid(
                "TSS context requires an explicit locus anchor strand",
            ));
        }
    };
    Ok((pair.0 as u64, pair.1 as u64))
}

fn genomic_strand(anchor: &GeneLocusGenomeAnchorBinding, local_reverse: bool) -> TssStrand {
    if local_reverse ^ (anchor.strand == Some('-')) {
        TssStrand::Minus
    } else {
        TssStrand::Plus
    }
}

fn span(
    geometry: &TssGeometry,
    start: u64,
    end: u64,
) -> Result<Option<TssContextSpan>, EngineError> {
    if start == 0 || end < start {
        return Err(invalid("Invalid genomic interval in TSS context"));
    }
    let left = start.max(geometry.start_1based);
    let right = end.min(geometry.end_1based);
    if right < left {
        return Ok(None);
    }
    let (local_start, local_end) = match geometry.strand {
        TssStrand::Plus => (
            left - geometry.start_1based,
            right - geometry.start_1based + 1,
        ),
        TssStrand::Minus => (geometry.end_1based - right, geometry.end_1based - left + 1),
    };
    Ok(Some(TssContextSpan {
        genomic_start_1based: start,
        genomic_end_1based: end,
        start_0based: local_start as usize,
        end_0based_exclusive: local_end as usize,
        clipped: left != start || right != end,
    }))
}

fn local_span(
    geometry: &TssGeometry,
    anchor: &GeneLocusGenomeAnchorBinding,
    start: usize,
    end: usize,
) -> Result<Option<TssContextSpan>, EngineError> {
    let (start, end) = local_genomic(anchor, start, end)?;
    span(geometry, start, end)
}

fn locus_sequence(bytes: &[u8], binding: &GeneLocusSequenceBinding) -> Result<String, EngineError> {
    let reader = bio::io::fasta::Reader::new(bytes);
    let mut records = reader.records();
    let record = records
        .next()
        .ok_or_else(|| invalid("Empty context locus FASTA"))?
        .map_err(|e| invalid(format!("Context FASTA: {e}")))?;
    record.check().map_err(|e| invalid(e.to_string()))?;
    if records.next().is_some() {
        return Err(invalid(
            "Context locus FASTA must contain exactly one full source sequence",
        ));
    }
    let sequence = std::str::from_utf8(record.seq()).map_err(|e| invalid(e.to_string()))?;
    if sequence.len() != binding.sequence_length_bp
        || !sequence.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
        || sha256_hex_bytes(sequence.as_bytes()) != digest(&binding.sequence_sha256)?
    {
        return Err(invalid(
            "Context FASTA does not match the locus report's exact loaded-sequence binding",
        ));
    }
    Ok(sequence.into())
}

fn verify_tata(
    report: &TataBoxScreenReport,
    locus: &GeneLocusEvidenceDisplayReport,
    anchor: &GeneLocusGenomeAnchorBinding,
) -> Result<(), EngineError> {
    let binding = locus
        .sequence_binding
        .as_ref()
        .ok_or_else(|| invalid("Unbound locus"))?;
    let mut unhashed = report.clone();
    unhashed.content_sha256.clear();
    let bytes = serde_json::to_vec(&unhashed).map_err(|e| invalid(e.to_string()))?;
    if report.schema != TATA_BOX_SCREEN_SCHEMA
        || report.request.seq_id != locus.seq_id
        || digest(&report.sequence_sha256)? != digest(&binding.sequence_sha256)?
        || report.genome_id.as_deref() != Some(anchor.genome_id.as_str())
        || report.chromosome.as_deref().map(chromosome) != Some(chromosome(&anchor.chromosome))
        || report.genomic_start_1based != Some(anchor.start_1based)
        || report.genomic_end_1based != Some(anchor.end_1based)
        || report.genomic_reverse != Some(anchor.strand == Some('-'))
        || digest(&report.content_sha256)? != sha256_hex_bytes(&bytes)
    {
        return Err(invalid(
            "TATA context schema, content digest or locus sequence/reference binding mismatch",
        ));
    }
    Ok(())
}

fn project(
    record: &TssRecord,
    reference: &TssReference,
    locus: &GeneLocusEvidenceDisplayReport,
    sequence: &str,
    locus_sha: &str,
    bindings: Vec<TssInputBinding>,
    tata: Option<(&TataBoxScreenReport, &str)>,
) -> Result<TssDetailContext, EngineError> {
    gentle_engine::tss_profiles::validate_record(record)?;
    let binding = locus.sequence_binding.as_ref().ok_or_else(|| {
        invalid("Legacy locus report has no sequence binding; regenerate it before TSS composition")
    })?;
    let anchor = binding
        .genome_anchor
        .as_ref()
        .ok_or_else(|| invalid("TSS context requires a genome-anchored locus"))?;
    let evidence = &locus.isoform_evidence;
    let g = &record.geometry;
    if anchor.start_1based == 0
        || anchor
            .end_1based
            .checked_sub(anchor.start_1based)
            .and_then(|n| n.checked_add(1))
            != Some(sequence.len())
        || !matches!(anchor.strand, Some('+') | Some('-'))
        || anchor.genome_id != reference.genome_id
        || ![reference.assembly.as_str(), reference.genome_id.as_str()]
            .contains(&evidence.assembly.as_str())
        || chromosome(&anchor.chromosome) != chromosome(&g.chromosome)
        || evidence.chromosome.as_deref().map(chromosome) != Some(chromosome(&g.chromosome))
        || locus.gene_symbol != record.gene_symbol
        || evidence.gene_symbol != record.gene_symbol
        || locus.gene_strand != g.strand.as_str()
        || reference
            .annotation_release
            .as_ref()
            .is_some_and(|r| evidence.annotation_release.as_ref() != Some(r))
        || g.start_1based < anchor.start_1based as u64
        || g.end_1based > anchor.end_1based as u64
    {
        return Err(invalid(format!(
            "TSS {}: locus gene, reference, release, strand or interval does not match",
            record.promoter_id
        )));
    }
    let (start, end) = if anchor.strand == Some('-') {
        (
            (anchor.end_1based as u64 - g.end_1based) as usize,
            (anchor.end_1based as u64 - g.start_1based + 1) as usize,
        )
    } else {
        (
            (g.start_1based - anchor.start_1based as u64) as usize,
            (g.end_1based - anchor.start_1based as u64 + 1) as usize,
        )
    };
    let local = sequence
        .get(start..end)
        .ok_or_else(|| invalid("Invalid TSS sequence slice"))?;
    let window_sequence = if (g.strand == TssStrand::Minus) != (anchor.strand == Some('-')) {
        GentleEngine::reverse_complement(local)
    } else {
        local.to_string()
    };
    if sha256_hex_bytes(window_sequence.as_bytes()) != record.sequence_sha256 {
        return Err(invalid(format!(
            "TSS {}: window sequence SHA-256 differs from the source locus slice",
            record.promoter_id
        )));
    }
    let mut context = TssDetailContext {
        schema: CONTEXT_SCHEMA.into(),
        promoter_id: record.promoter_id.clone(),
        geometry: g.clone(),
        window_sequence_sha256: record.sequence_sha256.clone(),
        locus_seq_id: locus.seq_id.clone(),
        locus_sequence_sha256: digest(&binding.sequence_sha256)?.into(),
        locus_report_sha256: locus_sha.into(),
        annotation_release: evidence.annotation_release.clone(),
        bindings,
        transcripts: vec![],
        occupancy: vec![],
        tata: None,
        warnings: locus.warnings.clone(),
        non_claims: CONTEXT_NON_CLAIMS.into(),
    };
    if let Some(splicing) = &evidence.splicing {
        let mut seen = BTreeSet::new();
        for transcript in &splicing.transcripts {
            if !record.transcripts.contains(&transcript.transcript_id) {
                continue;
            }
            if !seen.insert(transcript.transcript_id.clone()) {
                return Err(invalid("Ambiguous duplicate transcript in locus context"));
            }
            if !matches!(transcript.strand.as_str(), "+" | "-") {
                return Err(invalid("Transcript context has unknown local strand"));
            }
            let strand = genomic_strand(anchor, transcript.strand == "-");
            if strand != g.strand {
                return Err(invalid("TSS/transcript strand mismatch"));
            }
            let mut row = TssContextTranscript {
                transcript_id: transcript.transcript_id.clone(),
                label: transcript.label.clone(),
                genomic_strand: strand,
                exons: vec![],
                cds: vec![],
                codons: vec![],
            };
            let mut exons = transcript.exons.iter().collect::<Vec<_>>();
            exons.sort_by_key(|e| (e.start_1based, e.end_1based));
            if transcript.strand == "-" {
                exons.reverse();
            }
            for (i, exon) in exons.iter().enumerate() {
                if let Some(span) = local_span(g, anchor, exon.start_1based, exon.end_1based)? {
                    row.exons.push(TssContextExon {
                        number_5prime_to_3prime: i + 1,
                        span,
                    });
                }
            }
            let metrics = locus
                .transcript_metrics
                .iter()
                .filter(|m| m.transcript_id == transcript.transcript_id)
                .collect::<Vec<_>>();
            if metrics.len() > 1 {
                return Err(invalid("Ambiguous transcript CDS metrics"));
            }
            if let Some(metrics) = metrics.first() {
                if matches!(
                    metrics.coding_status.as_str(),
                    "complete_cds" | "partial_cds"
                ) {
                    for &(start, end) in &metrics.cds_ranges_local_1based {
                        if let Some(span) = local_span(g, anchor, start, end)? {
                            row.cds.push(span);
                        }
                    }
                } else if !metrics.cds_ranges_local_1based.is_empty() {
                    context.warnings.push(format!(
                        "{}: inferred or unclassified coding ranges are not displayed as annotated CDS ({})",
                        transcript.transcript_id, metrics.coding_status
                    ));
                }
            } else {
                context.warnings.push(format!(
                    "{}: CDS annotation unavailable",
                    transcript.transcript_id
                ));
            }
            for marker in locus
                .codon_markers
                .iter()
                .filter(|m| m.transcript_id == transcript.transcript_id)
            {
                let (position, _) = local_genomic(
                    anchor,
                    marker.local_position_1based,
                    marker.local_position_1based,
                )?;
                if position != marker.genomic_position_1based as u64
                    || (!marker.genomic_strand.is_empty()
                        && marker.genomic_strand != strand.as_str())
                    || marker.basis.trim().is_empty()
                {
                    return Err(invalid(
                        "Inconsistent annotated translation marker in locus context",
                    ));
                }
                if let Some(at) = span(g, position, position)? {
                    row.codons.push(TssContextCodon {
                        kind: marker.kind,
                        position_0based: at.start_0based,
                        genomic_position_1based: position,
                        basis: marker.basis.clone(),
                    });
                }
            }
            context.transcripts.push(row);
        }
        for id in &record.transcripts {
            if !seen.contains(id) {
                context.warnings.push(format!(
                    "{id}: transcript geometry unavailable in supplied locus report"
                ));
            }
        }
    } else {
        context
            .warnings
            .push("Transcript geometry unavailable in supplied locus report".into());
    }
    context
        .transcripts
        .sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));
    for group in &locus.occupancy_groups {
        for lane in &group.lanes {
            let mut row = TssContextOccupancyLane {
                group_id: group.group_id.clone(),
                group_label: group.label.clone(),
                scale_mode: group.scale_mode,
                lane_id: lane.lane.lane_id.clone(),
                label: lane
                    .display_label
                    .as_ref()
                    .unwrap_or(&lane.lane.display_label)
                    .clone(),
                state: lane.state,
                role: lane.role,
                source_id: lane.source_id.clone(),
                source_sha256: lane.source_sha256.clone(),
                source_kind: lane.lane.source_kind.clone(),
                condition: lane.condition_label.clone(),
                cell_line: lane.cell_line_label.clone(),
                assay: lane.assay.clone(),
                mark: lane.mark.clone(),
                factor: lane.factor.clone(),
                display_abs_max_score: lane.display_abs_max_score,
                intervals: vec![],
            };
            if lane.state == GeneLocusOccupancyLaneState::Available {
                if lane
                    .source_assembly
                    .as_ref()
                    .is_some_and(|a| a != &reference.assembly && a != &reference.genome_id)
                {
                    return Err(invalid(
                        "Available occupancy lane has an incompatible source assembly",
                    ));
                }
                for interval in &lane.lane.intervals {
                    let (start, end) = local_genomic(
                        anchor,
                        interval.local_start_1based,
                        interval.local_end_1based,
                    )?;
                    if start != interval.genomic_start_1based as u64
                        || end != interval.genomic_end_1based as u64
                    {
                        return Err(invalid(
                            "Occupancy genomic/local coordinates disagree with the source anchor",
                        ));
                    }
                    if let Some(span) = span(g, start, end)? {
                        row.intervals.push(TssContextSignal {
                            interval_id: interval.interval_id.clone(),
                            span,
                            score: interval.score,
                            label: interval.label.clone(),
                        });
                    }
                }
            }
            row.intervals.sort_by(|a, b| {
                (
                    a.span.start_0based,
                    a.span.end_0based_exclusive,
                    &a.interval_id,
                )
                    .cmp(&(
                        b.span.start_0based,
                        b.span.end_0based_exclusive,
                        &b.interval_id,
                    ))
            });
            context.occupancy.push(row);
        }
    }
    if let Some((report, hash)) = tata {
        verify_tata(report, locus, anchor)?;
        let mut rows = Vec::new();
        for evidence in &report.rows {
            if evidence.start_0based >= evidence.end_0based_exclusive {
                return Err(invalid("Empty TATA evidence interval"));
            }
            if let Some(span) = local_span(
                g,
                anchor,
                evidence.start_0based + 1,
                evidence.end_0based_exclusive,
            )? {
                rows.push(TssContextTataRow {
                    evidence: evidence.clone(),
                    span,
                    genomic_strand: genomic_strand(anchor, evidence.reverse),
                });
            }
        }
        rows.sort_by(|a, b| {
            (a.span.start_0based, &a.evidence.row_id)
                .cmp(&(b.span.start_0based, &b.evidence.row_id))
        });
        context.tata = Some(TssContextTata {
            report_id: report.report_id.clone(),
            report_sha256: hash.into(),
            motif_id: report.motif_id.clone(),
            matrix_sha256: report.matrix_sha256.clone(),
            score_policy: report.score_policy.clone(),
            epd_status: report.epd_status.clone(),
            rows,
            warnings: report.warnings.clone(),
        });
    }
    context.warnings.sort();
    context.warnings.dedup();
    Ok(context)
}

pub(super) fn attach(
    report: &mut TssProfileReport,
    path: &Path,
    should_continue: &mut dyn FnMut() -> bool,
) -> Result<(), EngineError> {
    let bytes = read(path, 2 * 1024 * 1024)?;
    let manifest: TssContextManifest =
        serde_json::from_slice(&bytes).map_err(|e| invalid(format!("Context manifest: {e}")))?;
    if manifest.schema != CONTEXT_INPUT_SCHEMA
        || manifest.reference != report.reference
        || manifest.genes.is_empty()
        || manifest.genes.len() > 256
    {
        return Err(invalid(
            "Context manifest requires the exact TSS reference and 1..256 explicit gene sources",
        ));
    }
    let base = path.parent().unwrap_or(Path::new("."));
    let mut sources = manifest.genes.iter().collect::<Vec<_>>();
    sources.sort_by(|a, b| a.gene_id.cmp(&b.gene_id));
    let mut seen = BTreeSet::new();
    let mut context_bytes = 0usize;
    for source in sources {
        if !should_continue() {
            return Err(GentleEngine::tfbs_cancelled_error("TSS context binding"));
        }
        if !seen.insert(&source.gene_id)
            || !report
                .windows
                .iter()
                .any(|w| w.record.gene_id == source.gene_id)
        {
            return Err(invalid(
                "Duplicate or unmatched gene_id in context manifest",
            ));
        }
        let locus_bytes = bound_file(base, &source.locus_report)?;
        let document = LocusDocument::from_json(&locus_bytes).map_err(invalid)?;
        let locus = document.locus();
        let binding = locus
            .sequence_binding
            .as_ref()
            .ok_or_else(|| invalid("Legacy locus report lacks sequence binding; regenerate it"))?;
        let sequence = locus_sequence(&bound_file(base, &source.locus_fasta)?, binding)?;
        let tata: Option<TataBoxScreenReport> = source
            .tata_report
            .as_ref()
            .map(|file| {
                serde_json::from_slice(&bound_file(base, file)?)
                    .map_err(|e| invalid(format!("TATA context report: {e}")))
            })
            .transpose()?;
        let locus_hash = digest(&source.locus_report.sha256)?;
        let tata_hash = source
            .tata_report
            .as_ref()
            .map(|f| digest(&f.sha256))
            .transpose()?;
        let mut bindings = vec![TssInputBinding {
            role: "tss_detail_context_manifest".into(),
            name: format!("context-{}.json", sha256_hex_bytes(&bytes)),
            sha256: sha256_hex_bytes(&bytes),
        }];
        for (role, input) in [
            ("locus_report", Some(&source.locus_report)),
            ("locus_fasta", Some(&source.locus_fasta)),
            ("tata_report", source.tata_report.as_ref()),
        ] {
            if let Some(input) = input {
                let hash = digest(&input.sha256)?;
                bindings.push(TssInputBinding {
                    role: format!("tss_detail_{role}"),
                    name: format!("{role}-{hash}"),
                    sha256: hash.into(),
                });
            }
        }
        for window in report
            .windows
            .iter_mut()
            .filter(|w| w.record.gene_id == source.gene_id)
        {
            if !should_continue() {
                return Err(GentleEngine::tfbs_cancelled_error("TSS context projection"));
            }
            if window.detail_context.is_some() {
                return Err(invalid(
                    "TSS context is already attached; enrich the original score report rather than replacing evidence",
                ));
            }
            let context = project(
                &window.record,
                &report.reference,
                locus,
                &sequence,
                locus_hash,
                bindings.clone(),
                tata.as_ref().zip(tata_hash),
            )?;
            let bytes = serde_json::to_vec(&context).map_err(|e| invalid(e.to_string()))?;
            context_bytes = context_bytes
                .checked_add(bytes.len())
                .filter(|total| *total <= 32 * 1024 * 1024)
                .ok_or_else(|| {
                    invalid("Projected TSS contexts exceed 32 MiB; split the input profile report")
                })?;
            window.detail_context = Some(context);
        }
    }
    Ok(())
}
