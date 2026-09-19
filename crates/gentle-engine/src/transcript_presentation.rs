//! Exact physical geometry and source-preserving annotation presentation.

use gentle_protocol::transcript_presentation::*;
use std::collections::{BTreeMap, BTreeSet};

/// Read-only local geometry for all annotation-comparison frontends.
pub mod projection;

fn hash(value: &impl serde::Serialize) -> Result<String, String> {
    let bytes = serde_json::to_vec(value).map_err(|e| e.to_string())?;
    Ok(ring::digest::digest(&ring::digest::SHA256, &bytes)
        .as_ref()
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect())
}

fn digest(raw: &str) -> bool {
    raw.len() == 64
        && raw
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

fn valid_interval(i: TranscriptInterval) -> bool {
    i.start_1based > 0 && i.start_1based <= i.end_1based && i.end_1based <= i64::MAX as u64
}

/// Compare full exon chains, independently of CDS/phase availability. Callers must
/// not use this structural projection as a cDNA identity or primer-specificity test.
pub fn compare_sources(
    report: &TranscriptStructurePresentation,
) -> Result<Vec<TranscriptSourceComparisonRow>, String> {
    validate(report)?;
    let supplied: BTreeSet<_> = report.sources.iter().map(|s| s.provider).collect();
    let mut chains = BTreeMap::<String, Vec<&TranscriptStructureGroup>>::new();
    for group in &report.structure_groups {
        chains
            .entry(group.exon_chain_id.clone())
            .or_default()
            .push(group);
    }
    chains
        .into_iter()
        .map(|(exon_chain_id, groups)| {
            let members: BTreeSet<_> = groups
                .iter()
                .flat_map(|g| g.member_record_ids.iter().cloned())
                .collect();
            let mut ensembl = BTreeSet::new();
            let mut refseq = BTreeSet::new();
            for record in report
                .records
                .iter()
                .filter(|r| members.contains(&r.record_id))
            {
                let source = report
                    .sources
                    .iter()
                    .find(|s| s.source_id == record.structure.source_id)
                    .ok_or("Unresolved comparison source")?;
                match source.provider {
                    TranscriptProvider::Ensembl => {
                        ensembl.insert(record.structure.transcript_id.clone());
                    }
                    TranscriptProvider::RefSeq => {
                        refseq.insert(record.structure.transcript_id.clone());
                    }
                }
            }
            let membership = match (ensembl.is_empty(), refseq.is_empty()) {
                (false, false) => TranscriptSourceMembership::Both,
                (false, true) if supplied.contains(&TranscriptProvider::RefSeq) => {
                    TranscriptSourceMembership::EnsemblOnly
                }
                (true, false) if supplied.contains(&TranscriptProvider::Ensembl) => {
                    TranscriptSourceMembership::RefSeqOnly
                }
                (false, true) => TranscriptSourceMembership::EnsemblOtherUnassessed,
                (true, false) => TranscriptSourceMembership::RefSeqOtherUnassessed,
                _ => return Err("Empty exon-chain comparison".into()),
            };
            Ok(TranscriptSourceComparisonRow {
                exon_chain_id,
                membership,
                ensembl_transcript_ids: ensembl.into_iter().collect(),
                refseq_transcript_ids: refseq.into_iter().collect(),
                member_record_ids: members.into_iter().collect(),
                structure_ids: groups.iter().map(|g| g.structure_id.clone()).collect(),
            })
        })
        .collect()
}

/// Build identities from full, uncropped geometry. Clipping is presentation only.
pub fn build(
    assembly: &str,
    chromosome: &str,
    locus_sequence_sha256: &str,
    mut sources: Vec<TranscriptSourceBinding>,
    mut records: Vec<SourceTranscriptStructure>,
) -> Result<TranscriptStructurePresentation, String> {
    if assembly.trim().is_empty()
        || chromosome.trim().is_empty()
        || !digest(locus_sequence_sha256)
        || sources.is_empty()
        || sources.len() > 16
        || records.len() > 256
    {
        return Err("Invalid structure presentation identity or source/record budget".into());
    }
    sources.sort_by(|a, b| a.source_id.cmp(&b.source_id));
    let mut source_ids = BTreeSet::new();
    for s in &sources {
        if s.source_id.is_empty()
            || !source_ids.insert(s.source_id.as_str())
            || s.assembly != assembly
            || s.locus_sequence_sha256 != locus_sequence_sha256
            || !digest(&s.annotation_sha256)
            || s.release.trim().is_empty()
            || s.accession.trim().is_empty()
            || s.chromosome.is_empty()
        {
            return Err(
                "Invalid, duplicated or assembly/sequence-mismatched annotation source binding"
                    .into(),
            );
        }
    }
    records.sort_by(|a, b| (&a.source_id, &a.transcript_id).cmp(&(&b.source_id, &b.transcript_id)));
    let mut record_keys = BTreeSet::new();
    let mut exons: BTreeMap<String, PhysicalExon> = BTreeMap::new();
    let mut groups: BTreeMap<String, TranscriptStructureGroup> = BTreeMap::new();
    let mut ticks: BTreeMap<(TranscriptProvider, u64, i8), Vec<String>> = BTreeMap::new();
    let mut rows = Vec::new();
    for mut r in records {
        let s = sources
            .iter()
            .find(|s| s.source_id == r.source_id)
            .ok_or("Unresolved transcript source")?;
        if r.transcript_id.trim().is_empty()
            || !record_keys.insert((r.source_id.clone(), r.transcript_id.clone()))
            || !matches!(r.strand, -1 | 1)
            || r.exons.is_empty()
            || r.exons.len() > 1024
        {
            return Err("Invalid or duplicate transcript identity/strand/exon geometry".into());
        }
        r.designations.sort();
        r.designations.dedup();
        r.notes.sort();
        r.notes.dedup();
        let chain: Vec<_> = r.exons.iter().map(|e| e.interval).collect();
        if chain.iter().any(|i| !valid_interval(*i))
            || chain.windows(2).any(|w| {
                if r.strand == 1 {
                    w[0].end_1based >= w[1].start_1based
                } else {
                    w[1].end_1based >= w[0].start_1based
                }
            })
        {
            return Err("Exon chain is overlapping, invalid or not in transcription order".into());
        }
        if let Some(cds) = &r.cds {
            if cds.len() > 1024
                || cds.iter().any(|c| {
                    !valid_interval(c.interval)
                        || c.phase.is_some_and(|p| p > 2)
                        || !chain.iter().any(|e| {
                            e.start_1based <= c.interval.start_1based
                                && c.interval.end_1based <= e.end_1based
                        })
                })
                || cds.windows(2).any(|w| {
                    if r.strand == 1 {
                        w[0].interval.end_1based >= w[1].interval.start_1based
                    } else {
                        w[1].interval.end_1based >= w[0].interval.start_1based
                    }
                })
            {
                return Err("Invalid, unordered or non-exonic CDS geometry".into());
            }
        }
        let record_id = format!(
            "transcript_record_{}",
            hash(&(&r.source_id, &r.transcript_id))?
        );
        let content_sha256 = hash(&r)?;
        let exon_ids = chain
            .iter()
            .map(|i| {
                let exon_id = format!(
                    "physical_exon_{}",
                    hash(&(assembly, chromosome, r.strand, i))?
                );
                exons
                    .entry(exon_id.clone())
                    .or_insert_with(|| PhysicalExon {
                        exon_id: exon_id.clone(),
                        interval: *i,
                        strand: r.strand,
                        member_record_ids: vec![],
                    })
                    .member_record_ids
                    .push(record_id.clone());
                Ok(exon_id)
            })
            .collect::<Result<Vec<_>, String>>()?;
        let exon_chain_id = format!(
            "exon_chain_{}",
            hash(&(assembly, chromosome, r.strand, &chain))?
        );
        let cds_geometry_id = format!(
            "cds_geometry_{}",
            hash(&(assembly, chromosome, r.strand, &r.cds))?
        );
        let structure_id = format!("structure_{}", hash(&(&exon_chain_id, &cds_geometry_id))?);
        groups
            .entry(structure_id.clone())
            .or_insert_with(|| TranscriptStructureGroup {
                structure_id,
                exon_chain_id,
                cds_geometry_id,
                exon_ids,
                cds: r.cds.clone(),
                member_record_ids: vec![],
            })
            .member_record_ids
            .push(record_id.clone());
        let first = chain[0];
        let start = if r.strand == 1 {
            first.start_1based
        } else {
            first.end_1based
        };
        ticks
            .entry((s.provider, start, r.strand))
            .or_default()
            .push(record_id.clone());
        rows.push(TranscriptPresentationRecord {
            record_id,
            content_sha256,
            structure: r,
        });
    }
    let tss_ticks: Vec<_> = ticks
        .iter()
        .map(|(&(provider, position, strand), members)| SourceTssTick {
            provider,
            genomic_position_1based: position,
            strand,
            member_record_ids: members.clone(),
            exact_cross_source_agreement: ticks
                .keys()
                .any(|&(p, g, s)| p != provider && g == position && s == strand),
        })
        .collect();
    let mut tss_deltas = Vec::new();
    for e in tss_ticks
        .iter()
        .filter(|t| t.provider == TranscriptProvider::Ensembl)
    {
        for r in tss_ticks
            .iter()
            .filter(|t| t.provider == TranscriptProvider::RefSeq && t.strand == e.strand)
        {
            tss_deltas.push(SourceTssDelta {
                ensembl_position_1based: e.genomic_position_1based,
                refseq_position_1based: r.genomic_position_1based,
                strand: e.strand,
                transcript_oriented_delta_bp: (r.genomic_position_1based as i64
                    - e.genomic_position_1based as i64)
                    * i64::from(e.strand),
            });
        }
    }
    let mut physical_exons: Vec<_> = exons.into_values().collect();
    physical_exons.sort_by_key(|e| (e.strand, e.interval));
    let mut result = TranscriptStructurePresentation {
        schema: SCHEMA.into(),
        content_sha256: String::new(),
        assembly: assembly.into(),
        chromosome: chromosome.into(),
        locus_sequence_sha256: locus_sequence_sha256.into(),
        sources,
        records: rows,
        physical_exons,
        structure_groups: groups.into_values().collect(),
        tss_ticks,
        tss_deltas,
        non_claims: NON_CLAIMS.into(),
    };
    result.content_sha256 = hash(&result)?;
    Ok(result)
}

/// Never trust serialized derived identities. Rebuild them from their bound records.
pub fn validate(report: &TranscriptStructurePresentation) -> Result<(), String> {
    let expected = build(
        &report.assembly,
        &report.chromosome,
        &report.locus_sequence_sha256,
        report.sources.clone(),
        report.records.iter().map(|r| r.structure.clone()).collect(),
    )?;
    if *report != expected {
        return Err("Transcript presentation content/geometry binding mismatch".into());
    }
    Ok(())
}

pub fn payload_coverage(requested: &[String], included: &[String]) -> TranscriptPayloadCoverage {
    let requested: BTreeSet<_> = requested.iter().cloned().collect();
    let included: BTreeSet<_> = included
        .iter()
        .filter(|id| requested.contains(*id))
        .cloned()
        .collect();
    let unassessed: Vec<_> = requested.difference(&included).cloned().collect();
    TranscriptPayloadCoverage {
        statement: format!(
            "{} linked transcript structures were not included in this context payload; status unassessed, not missing. {} of {} linked records included.",
            unassessed.len(),
            included.len(),
            requested.len()
        ),
        requested_transcript_ids: requested.into_iter().collect(),
        included_transcript_ids: included.into_iter().collect(),
        unassessed_transcript_ids: unassessed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Hand-crafted, synthetic genomic coordinates and accessions only. Recreated
    // by these functions; used for identity, strand and failure regression tests.
    fn source(id: &str, provider: TranscriptProvider) -> TranscriptSourceBinding {
        TranscriptSourceBinding {
            source_id: id.into(),
            provider,
            assembly: "synthetic-1".into(),
            release: "release-1".into(),
            accession: format!("resource-{id}"),
            chromosome: "test".into(),
            annotation_sha256: "a".repeat(64),
            locus_sequence_sha256: "b".repeat(64),
        }
    }
    fn record(source_id: &str, id: &str, strand: i8) -> SourceTranscriptStructure {
        let mut exons = vec![
            TranscriptExon {
                interval: TranscriptInterval {
                    start_1based: 100,
                    end_1based: 200,
                },
                source_exon_id: Some(format!("{id}-exon1")),
            },
            TranscriptExon {
                interval: TranscriptInterval {
                    start_1based: 300,
                    end_1based: 400,
                },
                source_exon_id: None,
            },
        ];
        if strand == -1 {
            exons.reverse();
        }
        SourceTranscriptStructure {
            source_id: source_id.into(),
            transcript_id: id.into(),
            label: id.into(),
            strand,
            exons,
            cds: Some(vec![TranscriptCds {
                interval: TranscriptInterval {
                    start_1based: 320,
                    end_1based: 370,
                },
                phase: Some(0),
            }]),
            designations: vec![],
            notes: vec![],
        }
    }
    fn report(records: Vec<SourceTranscriptStructure>) -> TranscriptStructurePresentation {
        build(
            "synthetic-1",
            "test",
            &"b".repeat(64),
            vec![
                source("e", TranscriptProvider::Ensembl),
                source("r", TranscriptProvider::RefSeq),
            ],
            records,
        )
        .unwrap()
    }

    fn projection_binding(strand: char) -> gentle_protocol::GeneLocusSequenceBinding {
        gentle_protocol::GeneLocusSequenceBinding {
            sequence_sha256: "b".repeat(64),
            sequence_length_bp: 500,
            genome_anchor: Some(gentle_protocol::GeneLocusGenomeAnchorBinding {
                genome_id: "synthetic-1".into(),
                chromosome: "test".into(),
                start_1based: 1,
                end_1based: 500,
                strand: Some(strand),
            }),
        }
    }

    #[test]
    fn local_projection_uses_anchor_orientation_once_and_retains_phase_provenance() {
        for gene_strand in [1, -1] {
            for anchor_strand in ['+', '-'] {
                let e = record("e", "E.1", gene_strand);
                let mut r = record("r", "R.4", gene_strand);
                r.cds.as_mut().unwrap()[0].phase = Some(1);
                let p = report(vec![e, r]);
                let local = projection::project(&p, &projection_binding(anchor_strand)).unwrap();
                assert_eq!(local.sources, p.sources);
                assert_eq!(local.rows.len(), 2);
                for row in &local.rows {
                    assert_eq!(row.comparison.membership, TranscriptSourceMembership::Both);
                    assert_eq!(
                        row.local_strand,
                        gene_strand * if anchor_strand == '-' { -1 } else { 1 }
                    );
                    let expected_start = match (gene_strand, anchor_strand) {
                        (1, '+') => 100,
                        (1, '-') => 401,
                        (-1, '+') => 400,
                        _ => 101,
                    };
                    assert_eq!(row.annotated_start_1based, Some(expected_start));
                    let cds = &row.cds.as_ref().unwrap()[0];
                    assert_eq!(
                        cds.interval,
                        if anchor_strand == '-' {
                            TranscriptInterval {
                                start_1based: 131,
                                end_1based: 181,
                            }
                        } else {
                            TranscriptInterval {
                                start_1based: 320,
                                end_1based: 370,
                            }
                        }
                    );
                    assert!(row.matches_local_chain(&row.exons, row.local_strand));
                    assert!(!row.matches_local_chain(&row.exons, -row.local_strand));
                    assert!(!row.matches_local_chain(&row.exons[..1], row.local_strand));
                    assert_eq!(
                        row.records[0].structure.cds.as_ref().unwrap()[0].phase,
                        cds.phase
                    );
                }
            }
        }
    }

    #[test]
    fn local_projection_clipping_never_creates_shared_or_loaded_chain_identity() {
        let e = record("e", "E", 1);
        let mut r = record("r", "R", 1);
        r.exons[0].interval.start_1based = 110;
        let p = report(vec![e, r]);
        let mut binding = projection_binding('+');
        binding.genome_anchor.as_mut().unwrap().start_1based = 150;
        binding.sequence_length_bp = 351;
        let local = projection::project(&p, &binding).unwrap();
        assert_eq!(local.rows[0].exons, local.rows[1].exons);
        assert_ne!(
            local.rows[0].comparison.exon_chain_id,
            local.rows[1].comparison.exon_chain_id
        );
        for row in &local.rows {
            assert_ne!(row.comparison.membership, TranscriptSourceMembership::Both);
            assert_eq!(row.annotated_start_1based, None);
            assert!(!row.matches_local_chain(&row.exons, row.local_strand));
        }
    }

    #[test]
    fn local_projection_rejects_stale_or_ambiguous_binding_and_keeps_unassessed() {
        let p = build(
            "synthetic-1",
            "test",
            &"b".repeat(64),
            vec![source("e", TranscriptProvider::Ensembl)],
            vec![record("e", "E", 1)],
        )
        .unwrap();
        let good = projection_binding('+');
        assert_eq!(
            projection::project(&p, &good).unwrap().rows[0]
                .comparison
                .membership,
            TranscriptSourceMembership::EnsemblOtherUnassessed
        );
        for which in 0..6 {
            let mut bad = good.clone();
            match which {
                0 => bad.sequence_sha256 = "a".repeat(64),
                1 => bad.sequence_length_bp += 1,
                2 => bad.genome_anchor.as_mut().unwrap().strand = None,
                3 => bad.genome_anchor.as_mut().unwrap().genome_id = "other".into(),
                4 => bad.genome_anchor.as_mut().unwrap().chromosome = "other".into(),
                _ => bad.genome_anchor.as_mut().unwrap().start_1based = 501,
            }
            assert!(projection::project(&p, &bad).is_err());
        }
    }

    #[test]
    fn identical_chains_merge_geometry_but_retain_every_source() {
        for strand in [1, -1] {
            let mut e = record("e", "E.1", strand);
            let mut r = record("r", "R.4", strand);
            e.designations.push(TranscriptDesignation {
                label: "Ensembl canonical".into(),
                field: "is_canonical".into(),
                value: "true".into(),
            });
            r.designations.push(TranscriptDesignation {
                label: "RefSeq Select".into(),
                field: "tag".into(),
                value: "RefSeq_Select".into(),
            });
            let p = report(vec![e.clone(), r.clone()]);
            assert_eq!(p.physical_exons.len(), 2);
            assert_eq!(p.structure_groups.len(), 1);
            assert_eq!(p.structure_groups[0].member_record_ids.len(), 2);
            assert_eq!(p.records.len(), 2);
            assert!(p.tss_ticks.iter().all(|t| t.exact_cross_source_agreement));
            assert_eq!(p.tss_deltas[0].transcript_oriented_delta_bp, 0);
            assert_eq!(
                serde_json::to_vec(&p).unwrap(),
                serde_json::to_vec(&report(vec![r, e])).unwrap()
            );
            validate(&p).unwrap();
        }
    }

    #[test]
    fn source_comparison_separates_full_chains_from_cds_and_missing_sources() {
        for strand in [1, -1] {
            let e = record("e", "E.1", strand);
            let mut r = record("r", "R.4", strand);
            r.cds.as_mut().unwrap()[0].phase = None;
            let p = report(vec![e.clone(), r.clone()]);
            let rows = compare_sources(&p).unwrap();
            assert_eq!(rows.len(), 1);
            assert_eq!(rows[0].membership, TranscriptSourceMembership::Both);
            assert_eq!(
                rows[0].structure_ids.len(),
                2,
                "unknown CDS phase is not identical CDS"
            );
            assert_eq!(rows[0].refseq_transcript_ids, vec!["R.4"]);
            r.exons[0].interval.start_1based += 1;
            let p = report(vec![e.clone(), r]);
            let rows = compare_sources(&p).unwrap();
            assert_eq!(
                rows.len(),
                2,
                "a shared exon or cropped overlap cannot establish chain identity"
            );
            assert!(
                rows.iter()
                    .any(|r| r.membership == TranscriptSourceMembership::EnsemblOnly)
            );
            assert!(
                rows.iter()
                    .any(|r| r.membership == TranscriptSourceMembership::RefSeqOnly)
            );
            let single = build(
                "synthetic-1",
                "test",
                &"b".repeat(64),
                vec![source("e", TranscriptProvider::Ensembl)],
                vec![e],
            )
            .unwrap();
            assert_eq!(
                compare_sources(&single).unwrap()[0].membership,
                TranscriptSourceMembership::EnsemblOtherUnassessed
            );
        }
    }

    #[test]
    fn shared_exon_does_not_merge_chains_or_different_cds() {
        let a = record("e", "A", 1);
        let mut b = record("r", "B", 1);
        let mut c = record("r", "C", 1);
        b.exons[0].interval.start_1based = 110;
        c.cds.as_mut().unwrap()[0].interval.start_1based = 330;
        let p = report(vec![a, b, c]);
        assert_eq!(p.physical_exons.len(), 3);
        assert_eq!(p.structure_groups.len(), 3);
        assert_eq!(
            p.structure_groups
                .iter()
                .map(|g| &g.exon_chain_id)
                .collect::<BTreeSet<_>>()
                .len(),
            2
        );
        assert_eq!(
            p.structure_groups
                .iter()
                .map(|g| &g.cds_geometry_id)
                .collect::<BTreeSet<_>>()
                .len(),
            2
        );
        assert_eq!(
            p.physical_exons
                .iter()
                .find(|e| e.interval.start_1based == 300)
                .unwrap()
                .member_record_ids
                .len(),
            3
        );
        assert!(p.tss_ticks.iter().any(|t| !t.exact_cross_source_agreement));
    }

    #[test]
    fn shifted_minus_tss_uses_transcript_oriented_delta() {
        let e = record("e", "E", -1);
        let mut r = record("r", "R", -1);
        r.exons[0].interval.end_1based = 410;
        let p = report(vec![e, r]);
        assert_eq!(p.tss_deltas[0].transcript_oriented_delta_bp, -10);
        assert!(!p.tss_ticks[0].exact_cross_source_agreement);
    }

    #[test]
    fn absent_cds_is_not_explicit_noncoding() {
        let mut e = record("e", "E", 1);
        e.cds = None;
        let mut r = record("r", "R", 1);
        r.cds = Some(vec![]);
        assert_eq!(report(vec![e, r]).structure_groups.len(), 2);
    }

    #[test]
    fn invalid_geometry_and_independent_bindings_fail_closed() {
        let p = report(vec![record("e", "E", 1), record("r", "R", 1)]);
        for case in 0..7 {
            let mut changed = p.clone();
            match case {
                0 => changed.sources[0].assembly = "wrong".into(),
                1 => changed.sources[0].locus_sequence_sha256 = "c".repeat(64),
                2 => changed.sources[0].annotation_sha256 = "d".repeat(64),
                3 => changed.records[0].structure.exons[0].interval.end_1based = 301,
                4 => changed.tss_deltas[0].transcript_oriented_delta_bp = 99,
                5 => changed.physical_exons[0].member_record_ids.clear(),
                _ => changed.sources[0].release = "other-release".into(),
            }
            assert!(validate(&changed).is_err(), "case {case}");
        }
    }

    #[test]
    fn omitted_payload_records_are_unassessed_not_warnings_or_absence() {
        let p = payload_coverage(&["a".into(), "b".into(), "c".into()], &["a".into()]);
        assert_eq!(p.unassessed_transcript_ids, vec!["b", "c"]);
        assert!(p.statement.starts_with("2 linked transcript structures"));
        assert!(p.statement.contains("unassessed, not missing"));
    }
}
