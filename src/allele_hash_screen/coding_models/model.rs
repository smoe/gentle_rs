//! Shared exon projection, GTF phase checks and table-1 translation for C1.

use super::*;
use crate::engine::PrimerPayloadTranscript;
use crate::observed_transcript_structures::annotation::VersionState;

fn accession(id: &str, versions: &BTreeSet<Option<String>>) -> CodingAccession {
    let suffix = id
        .rsplit_once('.')
        .filter(|(_, v)| !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()));
    let mut state = if versions.len() > 1 {
        VersionState::Inconsistent
    } else if let Some(Some(version)) = versions.first() {
        if version.is_empty() || !version.bytes().all(|b| b.is_ascii_digit()) {
            VersionState::Inconsistent
        } else {
            VersionState::Declared(version.clone())
        }
    } else {
        suffix.map_or(VersionState::Absent, |(_, v)| {
            VersionState::Declared(v.into())
        })
    };
    if let (Some((_, v)), VersionState::Declared(declared)) = (suffix, &state)
        && v != declared
    {
        state = VersionState::Inconsistent;
    }
    let versioned_accession = match &state {
        VersionState::Declared(v) => Some(format!("{}.{v}", suffix.map_or(id, |(base, _)| base))),
        _ => None,
    };
    CodingAccession {
        id: id.into(),
        version_state: state,
        versioned_accession,
    }
}

pub(super) fn proteins(
    raw: Option<&[u8]>,
    budgets: &CodingModelBudgets,
    meter: &mut Meter<'_>,
) -> Result<BTreeMap<String, String>, EngineError> {
    let mut out = BTreeMap::<String, String>::new();
    let Some(raw) = raw else {
        return Ok(out);
    };
    let mut current: Option<String> = None;
    for (n, line) in std::str::from_utf8(raw).map_err(fail)?.lines().enumerate() {
        meter.charge(line.len() + 1)?;
        if n >= budgets.max_rows || line.len() > budgets.max_line_bytes {
            return Err(fail("protein_budget"));
        }
        if let Some(header) = line.strip_prefix('>') {
            let id = header.split_whitespace().next().unwrap_or("");
            if !token(id) || out.insert(id.into(), String::new()).is_some() {
                return Err(fail("protein_fasta_identity"));
            }
            current = Some(id.into());
        } else if !line.is_empty() {
            if !line
                .bytes()
                .all(|c| b"ACDEFGHIKLMNPQRSTVWYBXZJUO*".contains(&c))
            {
                return Err(fail("protein_fasta_sequence"));
            }
            let sequence = out
                .get_mut(
                    current
                        .as_ref()
                        .ok_or_else(|| fail("protein_fasta_header"))?,
                )
                .unwrap();
            if sequence.len().saturating_add(line.len()) > budgets.max_transcript_bases {
                return Err(fail("protein_budget"));
            }
            sequence.push_str(line);
        }
    }
    if out.is_empty() {
        return Err(fail("protein_fasta_records"));
    }
    for sequence in out.values_mut() {
        if sequence.ends_with('*') {
            sequence.pop();
        }
        if sequence.is_empty() || sequence.contains('*') {
            return Err(fail("protein_fasta_stop"));
        }
    }
    Ok(out)
}

fn project(
    ranges: &[(usize, usize)],
    blocks: &[geometry::Block],
    strand: char,
    meter: &mut Meter<'_>,
) -> Result<Vec<(usize, usize)>, EngineError> {
    if ranges.windows(2).any(|w| w[0].1 > w[1].0) {
        return Err(fail("overlapping_segments"));
    }
    let mut out = Vec::new();
    for &(a, z) in ranges {
        if a >= z {
            return Err(fail("segment_geometry"));
        }
        let mut covered = 0usize;
        for block in blocks {
            meter.charge(1)?;
            if let Some(v) = geometry::project(block, (a, z), strand) {
                covered += v.1 - v.0;
                out.push(v);
            }
        }
        if covered != z - a {
            return Err(fail("segment_outside_exons"));
        }
    }
    out.sort_unstable();
    Ok(out)
}

fn interval(ranges: &[(usize, usize)]) -> Result<Option<(usize, usize)>, EngineError> {
    if ranges.windows(2).any(|w| w[0].1 != w[1].0) {
        return Err(fail("noncontiguous_cdna_feature"));
    }
    Ok(ranges
        .first()
        .zip(ranges.last())
        .map(|(first, last)| (first.0, last.1)))
}

pub(super) fn derive(
    request: &CodingTranscriptModelsRequest,
    locus: &str,
    annotation: &PrimerPayloadTranscript,
    rows: Rows,
    proteins: &BTreeMap<String, String>,
    meter: &mut Meter<'_>,
) -> Result<CodingTranscriptModel, EngineError> {
    let exons = &annotation.genomic_exons;
    let strand = rows.strand.ok_or_else(|| fail("strand"))?;
    let reference_end = request
        .reference
        .start_0based
        .checked_add(locus.len())
        .ok_or_else(|| fail("reference_interval"))?;
    if annotation.genomic_strand != Some(strand)
        || annotation.genomic_contig.as_deref() != Some(&request.reference.contig)
        || annotation.target_name != request.reference.contig
        || exons.is_empty()
        || exons.iter().any(|&(a, z)| {
            a >= z
                || a < request.reference.start_0based
                || z > reference_end
                || z > annotation.target_length
        })
        || exons.windows(2).any(|w| w[0].1 > w[1].0)
        || rows.exons.iter().copied().collect::<Vec<_>>() != *exons
        || rows
            .transcript_spans
            .iter()
            .any(|s| *s != (exons[0].0, exons[exons.len() - 1].1))
    {
        return Err(fail("exon_geometry_mismatch"));
    }
    let blocks = read_ends::exon_map(exons, strand)?;
    let length = blocks.last().unwrap().query.1;
    if length > request.budgets.max_transcript_bases {
        return Err(fail("transcript_base_budget"));
    }
    let mut cdna = Vec::with_capacity(length);
    for block in &blocks {
        meter.charge(block.query.1 - block.query.0)?;
        let slice = &locus.as_bytes()[block.target.0 - request.reference.start_0based
            ..block.target.1 - request.reference.start_0based];
        if strand == '+' {
            cdna.extend_from_slice(slice);
        } else {
            cdna.extend(
                slice
                    .iter()
                    .rev()
                    .map(|b| crate::iupac_code::IupacCode::letter_complement(*b)),
            );
        }
    }
    let mut output = CodingTranscriptModel {
        transcript: accession(&annotation.transcript_id, &rows.transcript_versions),
        gene_id: request.selector.gene_id.clone(),
        biotypes: rows.biotypes.iter().cloned().collect(),
        strand,
        genomic_exons_0based_half_open: exons.clone(),
        genomic_utr_segments: rows
            .utr
            .iter()
            .map(|(kind, spans)| (kind.clone(), spans.iter().copied().collect()))
            .collect(),
        transcript_completeness_tags: rows
            .tags
            .iter()
            .filter(|t| matches!(t.as_str(), "mRNA_start_NF" | "mRNA_end_NF"))
            .cloned()
            .collect(),
        tags: rows.tags.iter().cloned().collect(),
        coding: None,
        translation_status: "no_cds".into(),
        translation_reasons: Vec::new(),
    };
    for spans in output.genomic_utr_segments.values() {
        project(spans, &blocks, strand, meter)?;
    }
    if rows.cds.is_empty() {
        if !rows.start.is_empty() || !rows.stop.is_empty() {
            return Err(fail("codon_without_cds"));
        }
        output.translation_status = if rows
            .biotypes
            .iter()
            .any(|b| b.contains("coding") || b == "nonsense_mediated_decay")
        {
            "coding_not_defined"
        } else {
            "noncoding"
        }
        .into();
        return Ok(output);
    }
    let cds_genomic = rows.cds.keys().copied().collect::<Vec<_>>();
    let cds_projected = project(&cds_genomic, &blocks, strand, meter)?;
    let (start, mut end) = interval(&cds_projected)?.unwrap();
    let start_segments = rows.start.iter().copied().collect::<Vec<_>>();
    let stop_segments = rows.stop.iter().copied().collect::<Vec<_>>();
    let starts = interval(&project(&start_segments, &blocks, strand, meter)?)?;
    let stops = interval(&project(&stop_segments, &blocks, strand, meter)?)?;
    let mut reasons = Vec::<String>::new();
    if starts.is_none_or(|(a, z)| a != start || z - a != 3) || rows.tags.contains("cds_start_NF") {
        reasons.push("incomplete_start".into());
    }
    if let Some((a, z)) = starts
        && (a != start || z > end || z - a > 3)
    {
        return Err(fail("start_geometry"));
    }
    let stop_complete = stops.is_some_and(|(a, z)| z - a == 3);
    if let Some((a, z)) = stops {
        if z - a > 3 {
            return Err(fail("stop_geometry"));
        }
        if stop_complete {
            match request.provider.cds_row_convention {
                CdsRowConvention::StopCodonExcludedFromCdsRows if a == end => end = z,
                CdsRowConvention::StopCodonIncludedInCdsRows if z == end && a >= start => (),
                _ => return Err(fail("cds_row_convention_contradiction")),
            }
        } else {
            let fits = match request.provider.cds_row_convention {
                CdsRowConvention::StopCodonExcludedFromCdsRows => a == end,
                CdsRowConvention::StopCodonIncludedInCdsRows => {
                    a >= end.saturating_sub(3) && z <= end
                }
            };
            if !fits {
                return Err(fail("stop_geometry"));
            }
        }
    }
    if !stop_complete || rows.tags.contains("cds_end_NF") {
        reasons.push("incomplete_end".into());
    }
    let mut ordered = rows.cds.iter().collect::<Vec<_>>();
    if strand == '-' {
        ordered.reverse();
    }
    let initial_frame = ordered[0].1.unwrap_or(0) as usize;
    let mut consumed = (3 - initial_frame) % 3;
    let mut phase_valid = true;
    for (&(a, z), frame) in ordered {
        meter.charge(1)?;
        if *frame != Some(((3 - consumed % 3) % 3) as u8) {
            phase_valid = false;
        }
        consumed += z - a;
    }
    if initial_frame != 0 && !reasons.iter().any(|s| s == "incomplete_start") {
        reasons.push("incomplete_start".into());
    }
    let completeness = reasons.clone();
    if !phase_valid {
        reasons.push("phase_inconsistent".into());
    }
    let mut utr_segments = BTreeMap::new();
    for (kind, ranges) in &rows.utr {
        let ranges = ranges.iter().copied().collect::<Vec<_>>();
        let projected = project(&ranges, &blocks, strand, meter)?;
        if projected.iter().any(|&(a, z)| match kind.as_str() {
            "five_prime_utr" => z > start,
            "three_prime_utr" => a < end,
            _ => a < end && z > start,
        }) {
            return Err(fail("utr_geometry"));
        }
        utr_segments.insert(kind.clone(), ranges);
    }
    let coding = &cdna[start..end];
    if coding.len() % 3 != 0 {
        reasons.push("length_not_multiple_of_3".into());
    }
    if coding.iter().any(|b| !b"ACGT".contains(b)) {
        reasons.push("ambiguous_sequence".into());
    }
    let start_codon = (coding.len() >= 3).then(|| String::from_utf8(coding[..3].to_vec()).unwrap());
    let stop_codon = stop_complete.then(|| String::from_utf8(cdna[end - 3..end].to_vec()).unwrap());
    if start_codon.as_deref() != Some("ATG") {
        reasons.push("start_not_atg".into());
    }
    let mut translated = String::new();
    if !coding.iter().any(|b| !b"ACGT".contains(b)) {
        for (index, chunk) in coding.chunks_exact(3).enumerate() {
            meter.charge(1)?;
            let aa = crate::AMINO_ACIDS.codon2aa([chunk[0], chunk[1], chunk[2]], Some(1));
            if aa == crate::amino_acids::STOP_CODON {
                if (index + 1) * 3 != coding.len() {
                    reasons.push("internal_stop".into());
                }
            } else {
                translated.push(aa);
            }
        }
        if stop_complete && !matches!(stop_codon.as_deref(), Some("TAA" | "TAG" | "TGA")) {
            reasons.push("terminal_not_stop".into());
        }
    }
    if rows.protein_ids.len() > 1 {
        return Err(fail("protein_identity_conflict"));
    }
    let protein = rows
        .protein_ids
        .first()
        .map(|id| accession(id, &rows.protein_versions));
    let translation_verified = reasons.is_empty();
    let supplied = if request.protein_fasta.is_none() {
        "not_supplied"
    } else if let Some(p) = protein
        .as_ref()
        .and_then(|p| p.versioned_accession.as_ref())
        .and_then(|id| proteins.get(id))
    {
        if !translation_verified {
            "not_assessable"
        } else if p == &translated {
            "identical"
        } else {
            "differs"
        }
    } else {
        "missing"
    };
    if supplied == "differs" {
        reasons.push("supplied_protein_mismatch".into());
    }
    output.translation_status = reasons
        .first()
        .cloned()
        .unwrap_or_else(|| "complete_verified".into());
    let mut seen = BTreeSet::new();
    reasons.retain(|r| seen.insert(r.clone()));
    output.translation_reasons = reasons;
    let mut normalized_genomic_coding_segments = Vec::new();
    for block in &blocks {
        let (a, z) = (start.max(block.query.0), end.min(block.query.1));
        if a < z {
            normalized_genomic_coding_segments.push(if strand == '+' {
                (
                    block.target.0 + a - block.query.0,
                    block.target.0 + z - block.query.0,
                )
            } else {
                (
                    block.target.1 - (z - block.query.0),
                    block.target.1 - (a - block.query.0),
                )
            });
        }
    }
    normalized_genomic_coding_segments.sort_unstable();
    output.coding = Some(CodingGeometry {
        segments: rows
            .cds
            .into_iter()
            .map(|(genomic_0based_half_open, gtf_frame)| CodingSegment {
                genomic_0based_half_open,
                gtf_frame,
            })
            .collect(),
        start_codon_segments: start_segments,
        stop_codon_segments: stop_segments,
        utr_segments,
        normalized_coding_cdna_0based_half_open: (start, end),
        normalized_genomic_coding_segments,
        cdna_cds_start_1based: start + 1,
        cdna_cds_end_1based: end,
        start_codon,
        stop_codon,
        cds_completeness_reasons: completeness,
        protein,
        translated_protein_sha256: translation_verified.then(|| digest(translated.as_bytes())),
        supplied_protein_comparison: supplied.into(),
    });
    Ok(output)
}
