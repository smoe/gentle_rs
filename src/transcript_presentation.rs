//! Local, hash-checked annotation adapters. No sequence retrieval or preferred TSS.

use gentle_protocol::transcript_presentation::*;
use std::{collections::BTreeMap, io::Read, path::Path};

fn canonical_chr(raw: &str) -> &str {
    raw.strip_prefix("chr").unwrap_or(raw)
}
fn bare_hash(raw: &str) -> &str {
    raw.strip_prefix("sha256:").unwrap_or(raw)
}

/// Check source-coherent content at report ingestion, including its enclosing locus.
pub fn validate_locus(
    locus: &gentle_protocol::isoform_evidence::GeneLocusEvidenceDisplayReport,
) -> Result<(), String> {
    if let Some(p) = &locus.transcript_presentation {
        gentle_engine::transcript_presentation::validate(p)?;
        let binding = locus
            .sequence_binding
            .as_ref()
            .ok_or("Transcript presentation requires locus sequence binding")?;
        let anchor = binding
            .genome_anchor
            .as_ref()
            .ok_or("Transcript presentation requires a genomic anchor")?;
        if p.assembly != locus.isoform_evidence.assembly
            || canonical_chr(&p.chromosome) != canonical_chr(&anchor.chromosome)
            || p.locus_sequence_sha256 != bare_hash(&binding.sequence_sha256)
            || p.records.iter().any(|r| match locus.gene_strand.as_str() {
                "+" => r.structure.strand != 1,
                "-" => r.structure.strand != -1,
                _ => true,
            })
        {
            return Err("Transcript presentation does not bind the enclosing locus".into());
        }
    }
    Ok(())
}

/// Read and join explicitly supplied annotations onto one sequence binding.
pub fn load(
    inputs: &[TranscriptAnnotationSource],
    base: &Path,
    assembly: &str,
    chromosome: &str,
    sequence_sha256: &str,
) -> Result<Option<TranscriptStructurePresentation>, String> {
    if inputs.is_empty() {
        return Ok(None);
    }
    if inputs.len() > 16 {
        return Err("At most 16 transcript annotation sources may be joined".into());
    }
    let mut sources = Vec::new();
    let mut records = Vec::new();
    for input in inputs {
        if input.assembly != assembly
            || bare_hash(&input.locus_sequence_sha256) != bare_hash(sequence_sha256)
            || input.gene_ids.is_empty()
            || input.gene_ids.len() > 32
        {
            return Err(
                "Transcript annotation assembly, locus hash or explicit gene scope mismatch".into(),
            );
        }
        let file = std::fs::File::open(base.join(&input.path)).map_err(|e| e.to_string())?;
        if !file.metadata().map_err(|e| e.to_string())?.is_file() {
            return Err("Annotation must be a regular file".into());
        }
        let mut bytes = Vec::new();
        file.take(128 * 1024 * 1024 + 1)
            .read_to_end(&mut bytes)
            .map_err(|e| e.to_string())?;
        if bytes.len() > 128 * 1024 * 1024 {
            return Err("Annotation exceeds 128 MiB; supply a bounded gene annotation".into());
        }
        let sha = crate::digest_utils::sha256_hex_bytes(&bytes);
        if sha != bare_hash(&input.sha256) {
            return Err("Transcript annotation file SHA-256 mismatch".into());
        }
        let source_id = format!(
            "annotation_source_{}",
            crate::digest_utils::sha256_hex_bytes(
                &serde_json::to_vec(&(
                    input.provider,
                    &input.assembly,
                    &input.release,
                    &input.accession,
                    &input.chromosome,
                    &sha
                ))
                .map_err(|e| e.to_string())?
            )
        );
        let parsed = match input.format {
            TranscriptAnnotationFormat::Gff3 => gff3(&bytes, input, &source_id, chromosome)?,
            TranscriptAnnotationFormat::EnsemblGeneEntry => {
                ensembl(&bytes, input, &source_id, chromosome)?
            }
        };
        if parsed.is_empty() {
            return Err(
                "No transcript records resolved for the explicit annotation gene scope".into(),
            );
        }
        sources.push(TranscriptSourceBinding {
            source_id,
            provider: input.provider,
            assembly: input.assembly.clone(),
            release: input.release.clone(),
            accession: input.accession.clone(),
            chromosome: input.chromosome.clone(),
            annotation_sha256: sha,
            locus_sequence_sha256: bare_hash(sequence_sha256).into(),
        });
        records.extend(parsed);
    }
    gentle_engine::transcript_presentation::build(
        assembly,
        canonical_chr(chromosome),
        bare_hash(sequence_sha256),
        sources,
        records,
    )
    .map(Some)
}

fn order(r: &mut SourceTranscriptStructure) {
    r.exons.sort_by_key(|e| e.interval);
    if let Some(cds) = &mut r.cds {
        cds.sort_by_key(|c| c.interval);
        if r.strand == -1 {
            cds.reverse();
        }
    }
    if r.strand == -1 {
        r.exons.reverse();
    }
}

fn designation(field: &str, value: &str) -> Option<TranscriptDesignation> {
    let label = match value {
        "Ensembl_canonical" => "Ensembl canonical",
        "MANE_Select" | "MANE Select" => "MANE Select",
        "MANE_Plus_Clinical" | "MANE Plus Clinical" => "MANE Plus Clinical",
        "RefSeq_Select" | "RefSeq Select" => "RefSeq Select",
        _ => return None,
    };
    Some(TranscriptDesignation {
        label: label.into(),
        field: field.into(),
        value: value.into(),
    })
}

fn gff3(
    bytes: &[u8],
    input: &TranscriptAnnotationSource,
    source_id: &str,
    chromosome: &str,
) -> Result<Vec<SourceTranscriptStructure>, String> {
    use bio::{bio_types::strand::Strand, io::gff};
    let text = std::str::from_utf8(bytes).map_err(|e| e.to_string())?;
    // Embedded FASTA is not an annotation row. The complete file is still hashed.
    let annotation = text.split("##FASTA").next().ok_or("Empty GFF3")?;
    for line in annotation.lines() {
        if let Some(build) = line.strip_prefix("#!genome-build ") {
            if build.trim() != input.assembly {
                return Err("GFF3 genome-build contradicts the declared assembly".into());
            }
        }
    }
    let all = gff::Reader::new(annotation.as_bytes(), gff::GffType::GFF3)
        .records()
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| e.to_string())?;
    let attr =
        |r: &gff::Record, key: &str| r.attributes().get_vec(key).cloned().unwrap_or_default();
    let selected: Vec<_> = all
        .iter()
        .filter(|r| r.seqname() == input.chromosome)
        .collect();
    let chr_ok = canonical_chr(&input.chromosome) == canonical_chr(chromosome)
        || selected.iter().any(|r| {
            r.feature_type() == "region"
                && attr(r, "chromosome")
                    .iter()
                    .any(|c| canonical_chr(c) == canonical_chr(chromosome))
        });
    if !chr_ok {
        return Err(
            "GFF3 reference accession has no source-declared chromosome mapping to the locus"
                .into(),
        );
    }
    let exon_parents: std::collections::BTreeSet<_> = selected
        .iter()
        .filter(|r| matches!(r.feature_type(), "exon" | "CDS"))
        .flat_map(|r| attr(r, "Parent"))
        .collect();
    let mut transcript_spans = BTreeMap::new();
    let mut transcripts = BTreeMap::new();
    for r in &selected {
        if matches!(r.feature_type(), "exon" | "CDS" | "gene" | "region")
            || !attr(r, "Parent").iter().any(|p| input.gene_ids.contains(p))
        {
            continue;
        }
        let id = attr(r, "ID")
            .first()
            .cloned()
            .ok_or("Transcript feature lacks ID")?;
        if !r.feature_type().ends_with("RNA")
            && !r.feature_type().ends_with("transcript")
            && !exon_parents.contains(&id)
        {
            continue;
        }
        transcript_spans.insert(id.clone(), (*r.start(), *r.end()));
        let transcript_id = attr(r, "transcript_id")
            .first()
            .cloned()
            .unwrap_or_else(|| {
                if input.provider == TranscriptProvider::RefSeq {
                    attr(r, "Name")
                        .first()
                        .cloned()
                        .unwrap_or_else(|| id.clone())
                } else {
                    id.strip_prefix("transcript:").unwrap_or(&id).into()
                }
            });
        let version = attr(r, "version").first().cloned();
        let transcript_id = if !transcript_id.contains('.') {
            version
                .map(|v| format!("{transcript_id}.{v}"))
                .unwrap_or(transcript_id)
        } else {
            transcript_id
        };
        let strand = match r.strand() {
            Some(Strand::Forward) => 1,
            Some(Strand::Reverse) => -1,
            _ => return Err("Transcript strand unresolved".into()),
        };
        let designations = attr(r, "tag")
            .iter()
            .filter_map(|v| designation("tag", v))
            .collect();
        let mut notes = vec![format!(
            "GFF3 ID={id}; Parent={}; feature_type={}",
            attr(r, "Parent").join(","),
            r.feature_type()
        )];
        for key in ["partial", "start_range", "end_range", "exception"] {
            for value in attr(r, key) {
                notes.push(format!("{key}={value}"));
            }
        }
        let explicitly_noncoding = matches!(
            r.feature_type(),
            "ncRNA" | "lnc_RNA" | "rRNA" | "tRNA" | "snRNA" | "snoRNA"
        );
        if transcripts
            .insert(
                id,
                SourceTranscriptStructure {
                    source_id: source_id.into(),
                    transcript_id,
                    label: attr(r, "Name").first().cloned().unwrap_or_default(),
                    strand,
                    exons: vec![],
                    cds: explicitly_noncoding.then(Vec::new),
                    designations,
                    notes,
                },
            )
            .is_some()
        {
            return Err("Duplicate GFF3 transcript ID".into());
        }
    }
    for r in selected {
        if !matches!(r.feature_type(), "exon" | "CDS") {
            continue;
        }
        for parent in attr(r, "Parent") {
            if let Some(t) = transcripts.get_mut(&parent) {
                let strand = match r.strand() {
                    Some(Strand::Forward) => 1,
                    Some(Strand::Reverse) => -1,
                    _ => 0,
                };
                if strand != t.strand {
                    return Err("Exon/CDS strand contradicts transcript".into());
                }
                let interval = TranscriptInterval {
                    start_1based: *r.start(),
                    end_1based: *r.end(),
                };
                if r.feature_type() == "exon" {
                    t.exons.push(TranscriptExon {
                        interval,
                        source_exon_id: attr(r, "ID").first().cloned(),
                    });
                } else {
                    t.cds.get_or_insert_with(Vec::new).push(TranscriptCds {
                        interval,
                        phase: r
                            .phase()
                            .clone()
                            .try_into()
                            .map_err(|()| "Invalid GFF3 CDS phase")?,
                    });
                }
            }
        }
    }
    let mut result = Vec::new();
    for (id, mut r) in transcripts {
        let bounds = (
            r.exons.iter().map(|e| e.interval.start_1based).min(),
            r.exons.iter().map(|e| e.interval.end_1based).max(),
        );
        if let Some(&(start, end)) = transcript_spans.get(&id) {
            if bounds != (Some(start), Some(end)) {
                return Err("GFF3 transcript bounds disagree with supplied exon-chain endpoints; incomplete/invalid structure payload".into());
            }
        }
        order(&mut r);
        result.push(r);
    }
    Ok(result)
}

fn ensembl(
    bytes: &[u8],
    input: &TranscriptAnnotationSource,
    source_id: &str,
    chromosome: &str,
) -> Result<Vec<SourceTranscriptStructure>, String> {
    let entry: crate::ensembl_gene::EnsemblGeneEntry =
        serde_json::from_slice(bytes).map_err(|e| e.to_string())?;
    if input.provider != TranscriptProvider::Ensembl
        || entry.schema != "gentle.ensembl_gene_entry.v1"
        || entry.assembly_name.as_deref() != Some(input.assembly.as_str())
        || entry.seq_region_name.as_deref() != Some(input.chromosome.as_str())
        || canonical_chr(&input.chromosome) != canonical_chr(chromosome)
        || !input.gene_ids.contains(&entry.gene_id)
    {
        return Err("Ensembl entry/source assembly, chromosome or gene binding mismatch".into());
    }
    entry.transcripts.into_iter().map(|t| {
        let mut designations=Vec::new();
        if t.is_canonical==Some(true) {designations.push(TranscriptDesignation {label:"Ensembl canonical".into(),field:"is_canonical".into(),value:"true".into()});}
        let mut r=SourceTranscriptStructure {source_id:source_id.into(),
            transcript_id:if t.transcript_id.contains('.') {t.transcript_id.clone()} else {t.transcript_version.map(|v|format!("{}.{v}",t.transcript_id)).unwrap_or(t.transcript_id.clone())},
            label:t.display_name.unwrap_or_default(),strand:t.strand.ok_or("Ensembl transcript strand unavailable")?,
            exons:Vec::new(),cds:None,designations,notes:vec!["Ensembl entry has no release field; release is the explicitly bound request declaration. CDS phases are not supplied by this entry format.".into()]};
        for exon in t.exons {
            if exon.seq_region_name.as_deref().is_some_and(|c|c!=input.chromosome) || exon.strand.is_some_and(|s|s!=r.strand) {return Err("Ensembl exon is on another reference/strand".into());}
            r.exons.push(TranscriptExon {interval:TranscriptInterval {start_1based:exon.start_1based as u64,end_1based:exon.end_1based as u64},source_exon_id:Some(exon.exon_version.map(|v|format!("{}.{v}",exon.exon_id)).unwrap_or(exon.exon_id))});
        }
        if let Some(c)=t.translation {
            if let (Some(start),Some(end))=(c.genomic_start_1based,c.genomic_end_1based) {
                if start>end || start==0 || !r.exons.iter().any(|e|e.interval.start_1based<=start as u64 && e.interval.end_1based>=start as u64)
                    || !r.exons.iter().any(|e|e.interval.start_1based<=end as u64 && e.interval.end_1based>=end as u64) {return Err("Ensembl translation endpoints are not exonic".into());}
                r.cds=Some(r.exons.iter().filter_map(|e| {
                    let a=e.interval.start_1based.max(start as u64);let b=e.interval.end_1based.min(end as u64);
                    (a<=b).then_some(TranscriptCds {interval:TranscriptInterval {start_1based:a,end_1based:b},phase:None})
                }).collect());
            }
        }
        order(&mut r); Ok(r)
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    // Entirely synthetic GFF3, regenerated in a temporary directory by each test.
    // Exercises the real parser and bound source loader, never private gene data.
    fn gff(strand: &str) -> String {
        format!(
            "##gff-version 3\n#!genome-build synthetic\nref.1\tRefSeq\tregion\t1\t900\t.\t+\t.\tID=ref.1;chromosome=test\nref.1\tRefSeq\tmRNA\t100\t400\t.\t{strand}\t.\tID=rna-A;Parent=gene-G;Name=SYN.1;tag=RefSeq_Select,MANE_Select\nref.1\tRefSeq\texon\t100\t200\t.\t{strand}\t.\tID=exon-1;Parent=rna-A\nref.1\tRefSeq\texon\t300\t400\t.\t{strand}\t.\tID=exon-2;Parent=rna-A\nref.1\tRefSeq\tCDS\t320\t370\t.\t{strand}\t0\tParent=rna-A\n"
        )
    }
    fn input(path: &Path, content: &str) -> TranscriptAnnotationSource {
        TranscriptAnnotationSource {
            path: path.to_string_lossy().into(),
            sha256: crate::digest_utils::sha256_hex_bytes(content.as_bytes()),
            provider: TranscriptProvider::RefSeq,
            format: TranscriptAnnotationFormat::Gff3,
            assembly: "synthetic".into(),
            release: "source-release".into(),
            accession: "source-accession".into(),
            chromosome: "ref.1".into(),
            locus_sequence_sha256: "a".repeat(64),
            gene_ids: vec!["gene-G".into()],
        }
    }
    #[test]
    fn transcript_source_gff3_preserves_designations_phase_and_source_reference_mapping() {
        for strand in ["+", "-"] {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("synthetic.gff3");
            let data = gff(strand);
            std::fs::write(&path, &data).unwrap();
            let p = load(
                &[input(&path, &data)],
                dir.path(),
                "synthetic",
                "test",
                &"a".repeat(64),
            )
            .unwrap()
            .unwrap();
            assert_eq!(p.records[0].structure.designations.len(), 2);
            assert_eq!(p.records[0].structure.transcript_id, "SYN.1");
            assert_eq!(p.sources[0].chromosome, "ref.1");
            assert_eq!(p.chromosome, "test");
            assert_eq!(
                p.records[0].structure.cds.as_ref().unwrap()[0].phase,
                Some(0)
            );
            assert_eq!(
                p.tss_ticks[0].genomic_position_1based,
                if strand == "+" { 100 } else { 400 }
            );
            gentle_engine::transcript_presentation::validate(&p).unwrap();
        }
    }
    #[test]
    fn transcript_source_independent_file_assembly_locus_and_gene_bindings_reject_tamper() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("synthetic.gff3");
        let data = gff("+");
        std::fs::write(&path, &data).unwrap();
        let original = input(&path, &data);
        for case in 0..6 {
            let mut changed = original.clone();
            match case {
                0 => changed.sha256 = "b".repeat(64),
                1 => changed.assembly = "other".into(),
                2 => changed.locus_sequence_sha256 = "b".repeat(64),
                3 => changed.gene_ids = vec!["unresolved".into()],
                4 => changed.chromosome = "wrong".into(),
                _ => changed.release.clear(),
            }
            assert!(
                load(&[changed], dir.path(), "synthetic", "test", &"a".repeat(64)).is_err(),
                "case {case}"
            );
        }
        let bad = data.replace("100\t200", "200\t100");
        std::fs::write(&path, &bad).unwrap();
        assert!(
            load(
                &[input(&path, &bad)],
                dir.path(),
                "synthetic",
                "test",
                &"a".repeat(64)
            )
            .is_err()
        );
    }
    #[test]
    fn transcript_source_ensembl_entry_geometry_does_not_invent_mane_or_cds_phases() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("entry.json");
        let data=serde_json::json!({"schema":"gentle.ensembl_gene_entry.v1","gene_id":"gene-G","assembly_name":"synthetic","seq_region_name":"test",
            "transcripts":[{"transcript_id":"ENSYN","transcript_version":3,"strand":1,"is_canonical":true,
                "exons":[{"exon_id":"ex1","start_1based":100,"end_1based":200,"seq_region_name":"test"},{"exon_id":"ex2","start_1based":300,"end_1based":400,"seq_region_name":"test"}],
                "translation":{"genomic_start_1based":180,"genomic_end_1based":370}}]}).to_string();
        std::fs::write(&path, &data).unwrap();
        let mut source = input(&path, &data);
        source.format = TranscriptAnnotationFormat::EnsemblGeneEntry;
        source.provider = TranscriptProvider::Ensembl;
        source.chromosome = "test".into();
        let p = load(&[source], dir.path(), "synthetic", "test", &"a".repeat(64))
            .unwrap()
            .unwrap();
        let r = &p.records[0].structure;
        assert_eq!(r.transcript_id, "ENSYN.3");
        assert_eq!(r.designations.len(), 1);
        assert_eq!(r.designations[0].label, "Ensembl canonical");
        assert_eq!(r.cds.as_ref().unwrap().len(), 2);
        assert!(r.cds.as_ref().unwrap().iter().all(|c| c.phase.is_none()));
    }

    #[test]
    fn transcript_source_gff3_ignores_gene_regulatory_children_but_rejects_truncated_chains() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("gene.gff3");
        let data = format!(
            "{}ref.1\tRefSeq\tregulatory_region\t80\t95\t.\t+\t.\tID=promoter;Parent=gene-G\n",
            gff("+")
        );
        std::fs::write(&path, &data).unwrap();
        let p = load(
            &[input(&path, &data)],
            dir.path(),
            "synthetic",
            "test",
            &"a".repeat(64),
        )
        .unwrap()
        .unwrap();
        assert_eq!(p.records.len(), 1);
        let truncated = data.replace("exon\t100\t200", "exon\t101\t200");
        std::fs::write(&path, &truncated).unwrap();
        assert!(
            load(
                &[input(&path, &truncated)],
                dir.path(),
                "synthetic",
                "test",
                &"a".repeat(64)
            )
            .unwrap_err()
            .contains("transcript bounds")
        );
    }
}
