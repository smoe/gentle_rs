//! Exact annotation identities and streamed FASTA verification with bounded edge contexts.

use super::*;

pub(super) struct Catalog {
    pub annotation: PrimerPayloadAnnotation,
    pub by_target: BTreeMap<String, Vec<usize>>,
    pub lengths: BTreeMap<String, usize>,
}

fn names(ids: &[String]) -> bool {
    !ids.is_empty()
        && ids
            .iter()
            .all(|s| !s.is_empty() && s.len() <= 512 && !s.chars().any(char::is_control))
        && ids.windows(2).all(|w| w[0] < w[1])
}

pub(super) fn catalog(
    reference: &PrimerPayloadReference,
    b: &PrimerPayloadBudgets,
) -> Result<Catalog, EngineError> {
    let annotation: PrimerPayloadAnnotation = read_json(&reference.annotation, b)?;
    if annotation.schema != PRIMER_PAYLOAD_ANNOTATION_SCHEMA
        || annotation.reference_id != reference.reference_id
        || annotation.assembly != reference.assembly
        || annotation.release != reference.release
        || annotation.reference_fasta_sha256 != reference.fasta.sha256
        || annotation.coordinate_convention != "zero_based_half_open"
        || annotation.transcripts.is_empty()
        || annotation.transcripts.len() > b.max_annotation_rows
    {
        return Err(bad("primer_payload_annotation_binding"));
    }
    let mut ids = BTreeSet::new();
    let mut by_target = BTreeMap::<String, Vec<usize>>::new();
    let mut lengths = BTreeMap::new();
    for (i, transcript) in annotation.transcripts.iter().enumerate() {
        if !names(&transcript.gene_ids)
            || transcript.transcript_id.is_empty()
            || !ids.insert(&transcript.transcript_id)
            || transcript.target_name.is_empty()
            || transcript.target_name.chars().any(char::is_whitespace)
            || transcript.target_length == 0
            || lengths
                .insert(transcript.target_name.clone(), transcript.target_length)
                .is_some_and(|n| n != transcript.target_length)
            || transcript.genomic_exons.windows(2).any(|w| w[0].1 > w[1].0)
            || transcript.genomic_exons.iter().any(|(a, z)| a >= z)
        {
            return Err(bad("primer_payload_annotation_identity_or_geometry"));
        }
        match reference.scope {
            PrimerPayloadReferenceScope::Genome => {
                if transcript.genomic_contig.as_deref() != Some(&transcript.target_name)
                    || !matches!(transcript.genomic_strand, Some('+' | '-'))
                    || transcript.genomic_exons.is_empty()
                    || transcript
                        .genomic_exons
                        .iter()
                        .any(|(_, z)| *z > transcript.target_length)
                {
                    return Err(bad("primer_payload_genomic_annotation_required"));
                }
            }
            PrimerPayloadReferenceScope::Transcriptome => {
                if transcript.genomic_exons.is_empty() {
                    if transcript.genomic_contig.is_some() || transcript.genomic_strand.is_some() {
                        return Err(bad("primer_payload_annotation_projection"));
                    }
                } else if transcript
                    .genomic_contig
                    .as_deref()
                    .is_none_or(str::is_empty)
                    || !matches!(transcript.genomic_strand, Some('+' | '-'))
                    || transcript
                        .genomic_exons
                        .iter()
                        .try_fold(0usize, |n, (a, z)| n.checked_add(z - a))
                        != Some(transcript.target_length)
                {
                    return Err(bad("primer_payload_annotation_projection"));
                }
            }
        }
        by_target
            .entry(transcript.target_name.clone())
            .or_default()
            .push(i);
    }
    Ok(Catalog {
        annotation,
        by_target,
        lengths,
    })
}

/// Clip reference/exon overlap back to original query coordinates, including reverse mappings.
pub(super) fn project(
    block: &paf::Block,
    range: (usize, usize),
    strand: char,
) -> Option<(usize, usize)> {
    let (a, z) = (block.target.0.max(range.0), block.target.1.min(range.1));
    if a >= z {
        return None;
    }
    Some(if strand == '+' {
        (
            block.query.0 + a - block.target.0,
            block.query.0 + z - block.target.0,
        )
    } else {
        (
            block.query.1 - (z - block.target.0),
            block.query.1 - (a - block.target.0),
        )
    })
}

pub(super) fn annotate(
    alignment: &mut paf::Alignment,
    reference: &PrimerPayloadReference,
    cat: &Catalog,
    work: &mut usize,
    b: &PrimerPayloadBudgets,
) -> Result<(), EngineError> {
    let seg = &mut alignment.segment;
    let mut genes = BTreeSet::new();
    let mut transcripts = BTreeSet::new();
    let mut splice = BTreeMap::<(String, usize, usize, char), BTreeSet<String>>::new();
    let indices = cat.by_target.get(&seg.target_name);
    if indices.is_none() && reference.scope == PrimerPayloadReferenceScope::Transcriptome {
        return Err(bad("primer_payload_unannotated_transcript"));
    }
    for &index in indices.into_iter().flatten() {
        let tx = &cat.annotation.transcripts[index];
        let comparisons = tx
            .genomic_exons
            .len()
            .max(1)
            .checked_mul(alignment.blocks.len().max(1) + alignment.skips.len())
            .ok_or_else(|| bad("primer_payload_annotation_comparison_budget"))?;
        *work = work
            .checked_add(comparisons)
            .filter(|n| *n <= b.max_chain_comparisons)
            .ok_or_else(|| bad("primer_payload_annotation_comparison_budget"))?;
        let overlaps = reference.scope == PrimerPayloadReferenceScope::Transcriptome
            || tx
                .genomic_exons
                .iter()
                .any(|(a, z)| *a < seg.target_interval.1 && *z > seg.target_interval.0);
        if !overlaps {
            continue;
        }
        genes.extend(tx.gene_ids.iter().cloned());
        transcripts.insert(tx.transcript_id.clone());
        if reference.scope == PrimerPayloadReferenceScope::Genome {
            for &(a, z) in &alignment.skips {
                if tx
                    .genomic_exons
                    .windows(2)
                    .any(|w| w[0].1 == a && w[1].0 == z)
                {
                    let strand = tx.genomic_strand.expect("validated");
                    let (donor, acceptor) = if strand == '+' { (a, z) } else { (z, a) };
                    splice
                        .entry((seg.target_name.clone(), donor, acceptor, strand))
                        .or_default()
                        .insert(tx.transcript_id.clone());
                }
            }
        }
    }
    seg.gene_ids = genes.into_iter().collect();
    seg.transcript_ids = transcripts.into_iter().collect();
    seg.splice_edges = splice
        .into_iter()
        .map(
            |((contig, donor, acceptor, strand), tx)| PrimerPayloadSpliceEdge {
                contig,
                donor,
                acceptor,
                strand,
                transcript_ids: tx.into_iter().collect(),
            },
        )
        .collect();
    Ok(())
}

pub(super) fn expected_spans(
    a: &paf::Alignment,
    expected: &[&PrimerPayloadExpectedTarget],
    cat: &Catalog,
    scope: PrimerPayloadReferenceScope,
    masks: &[(usize, usize)],
) -> Vec<(usize, usize)> {
    let mut spans = Vec::new();
    for &index in cat
        .by_target
        .get(&a.segment.target_name)
        .into_iter()
        .flatten()
    {
        let tx = &cat.annotation.transcripts[index];
        if !expected.iter().any(|e| {
            e.transcript_ids.contains(&tx.transcript_id)
                || (!e.gene_ids.is_empty() && tx.gene_ids.iter().any(|g| e.gene_ids.contains(g)))
        }) {
            continue;
        }
        for block in &a.blocks {
            if scope == PrimerPayloadReferenceScope::Transcriptome {
                spans.push(block.query);
            } else {
                for exon in &tx.genomic_exons {
                    if let Some(q) = project(block, *exon, a.segment.strand) {
                        spans.push(q);
                    }
                }
            }
        }
    }
    paf::subtract_masks(spans, masks)
}

/// A context key is reference target plus the actual half-open range; absent flanks stay unavailable.
pub(super) type ContextKey = (String, usize, usize);

pub(super) fn verify_fasta(
    reference: &PrimerPayloadReference,
    lengths: &BTreeMap<String, usize>,
    wanted: &BTreeSet<ContextKey>,
    b: &PrimerPayloadBudgets,
    ctl: &mut Control<'_>,
) -> Result<BTreeMap<ContextKey, Vec<u8>>, EngineError> {
    let total = wanted
        .iter()
        .try_fold(0usize, |n, (_, a, z)| n.checked_add(z - a))
        .filter(|n| *n <= b.max_reference_context_bases)
        .ok_or_else(|| bad("primer_payload_reference_context_budget"))?;
    let _ = total;
    let mut contexts = wanted
        .iter()
        .cloned()
        .map(|k| (k, Vec::new()))
        .collect::<BTreeMap<_, _>>();
    let mut ranges = BTreeMap::<String, Vec<(usize, usize)>>::new();
    for (target, a, z) in wanted {
        ranges.entry(target.clone()).or_default().push((*a, *z));
    }
    for values in ranges.values_mut() {
        values.sort_unstable();
    }
    let mut reader = BufReader::new(VerifiedReader::open(&reference.fasta, b.max_input_bytes)?);
    let (mut current, mut offset) = (String::new(), 0usize);
    let mut seen = BTreeSet::new();
    let mut active = Vec::new();
    let mut next_range = 0usize;
    let mut first_line = true;
    let finish = |name: &str, length: usize| -> Result<(), EngineError> {
        if lengths.get(name).is_some_and(|n| *n != length) {
            return Err(bad("primer_payload_fasta_target_length"));
        }
        Ok(())
    };
    while let Some(text) = line(&mut reader, b.max_line_bytes)? {
        if first_line && ctl.stop(true) {
            return Err(bad(INTERRUPTED));
        }
        first_line = false;
        ctl.check()?;
        if let Some(header) = text.strip_prefix('>') {
            if !current.is_empty() {
                finish(&current, offset)?;
            }
            current = header.split_ascii_whitespace().next().unwrap_or("").into();
            if current.is_empty()
                || !seen.insert(current.clone())
                || seen.len() > b.max_annotation_rows
            {
                return Err(bad("primer_payload_fasta_identity_or_budget"));
            }
            offset = 0;
            active.clear();
            next_range = 0;
        } else {
            if current.is_empty()
                || text.is_empty()
                || !text
                    .bytes()
                    .all(|x| b"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(&x))
            {
                return Err(bad("primer_payload_fasta_sequence"));
            }
            let end = offset
                .checked_add(text.len())
                .ok_or_else(|| bad("primer_payload_fasta_length_overflow"))?;
            if let Some(rs) = ranges.get(&current) {
                while next_range < rs.len() && rs[next_range].0 < end {
                    active.push(rs[next_range]);
                    next_range += 1;
                }
                active.retain(|(_, z)| *z > offset);
                for &(a, z) in &active {
                    let s = a.max(offset);
                    let e = z.min(end);
                    if s < e {
                        contexts
                            .get_mut(&(current.clone(), a, z))
                            .expect("admitted range")
                            .extend(
                                text.as_bytes()[s - offset..e - offset]
                                    .iter()
                                    .map(u8::to_ascii_uppercase),
                            );
                    }
                }
            }
            offset = end;
        }
    }
    if !current.is_empty() {
        finish(&current, offset)?;
    }
    reader.into_inner().finish()?;
    if !lengths.keys().all(|name| seen.contains(name))
        || contexts.iter().any(|((_, a, z), s)| s.len() != z - a)
    {
        return Err(bad("primer_payload_fasta_missing_target_or_context"));
    }
    Ok(contexts)
}
