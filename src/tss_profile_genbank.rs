//! Annotated TSS sequence projection from verified context and stored scores.
//! No sequence retrieval, feature inference, thresholding or motif rescoring.

use super::*;
use gb_io::seq::{After, Before, Feature, Location, Seq};

fn clean(text: &str) -> String {
    text.chars()
        .map(|c| if c.is_control() { ' ' } else { c })
        .collect()
}

fn feature(kind: &'static str, location: Location, label: &str, note: String) -> Feature {
    Feature {
        kind: kind.into(),
        location,
        qualifiers: vec![
            ("label".into(), Some(clean(label))),
            ("note".into(), Some(clean(&note))),
        ],
    }
}

fn location(span: &TssContextSpan, g: &TssGeometry) -> Location {
    let (before, after) = if g.strand == TssStrand::Plus {
        (
            span.genomic_start_1based < g.start_1based,
            span.genomic_end_1based > g.end_1based,
        )
    } else {
        (
            span.genomic_end_1based > g.end_1based,
            span.genomic_start_1based < g.start_1based,
        )
    };
    Location::Range(
        (span.start_0based as i64, Before(before)),
        (span.end_0based_exclusive as i64, After(after)),
    )
}

fn span_note(span: &TssContextSpan) -> String {
    format!(
        "Source genomic interval {}..{} (1-based inclusive); cropped_to_window={}",
        span.genomic_start_1based, span.genomic_end_1based, span.clipped
    )
}

pub(super) fn bytes(
    report: &TssProfileReport,
    window: &TssProfileWindow,
) -> Result<Vec<u8>, EngineError> {
    let context = window.detail_context.as_ref().ok_or_else(|| invalid("GenBank export requires verified TSS context; supply --context-manifest for every gene"))?;
    let bases = context.window_sequence.as_deref().ok_or_else(|| invalid("Legacy context has no stored bases; enrich the original score-only report with --context-manifest before GenBank export"))?;
    let r = &window.record;
    let g = &r.geometry;
    if Some(bases.len()) != g.length()
        || !bases.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
        || sha256_hex_bytes(bases.as_bytes()) != r.sequence_sha256
        || context.window_sequence_sha256 != r.sequence_sha256
        || context.geometry != *g
    {
        return Err(invalid(
            "GenBank bases or context do not match the exact TSS sequence binding",
        ));
    }
    let mut seq = Seq::empty();
    seq.name = Some(format!(
        "TSS_{}",
        &sha256_hex_bytes(r.promoter_id.as_bytes())[..12]
    ));
    seq.seq = bases.as_bytes().to_vec();
    seq.molecule_type = Some("DNA".into());
    seq.definition = Some(clean(&format!(
        "{} annotated TSS window {}; transcript-oriented genomic DNA (not a spliced transcript or reporter construct)",
        r.gene_symbol, r.promoter_id
    )));
    seq.comments = vec![
        "The LOCUS date 01-JAN-1970 is the writer's deterministic missing-date placeholder, not an annotation, experiment or sample date.".into(),
        format!("GENtle promoter_id={}; selected={}; gene_id={}; sequence_sha256={}", clean(&r.promoter_id), window.selected, clean(&r.gene_id), r.sequence_sha256),
        format!("Reference={}; assembly={}; annotation_release={}; chromosome={}; genomic={}..{}; genomic_strand={}; local_axis=transcript_5prime_to_3prime; TSS_local_1based={}", clean(&report.reference.genome_id), clean(&report.reference.assembly), report.reference.annotation_release.as_deref().unwrap_or("not supplied"), clean(&g.chromosome), g.start_1based, g.end_1based, g.strand.as_str(), g.upstream_bp + 1),
        format!("locus_report_sha256={}; locus_sequence_sha256={}", context.locus_report_sha256, context.locus_sequence_sha256),
        "Locations are local 1-based inclusive; < and > mark clipping. Original genomic spans remain in notes. CDS fragments are misc_features because this context does not bind coding phase. Signal intervals are not individual reads. Stored motif peaks are predictions, not all possible sites or measured binding.".into(),
        clean(&report.non_claims),
    ];
    seq.comments
        .extend(context.warnings.iter().map(|w| clean(w)));
    seq.features.push(feature(
        "source",
        Location::simple_range(0, bases.len() as i64),
        &r.gene_symbol,
        format!(
            "Gene {}; chromosome {}; genomic {}..{}; genomic strand {}; transcripts={}",
            clean(&r.gene_id),
            clean(&g.chromosome),
            g.start_1based,
            g.end_1based,
            g.strand.as_str(),
            r.transcripts
                .iter()
                .map(|s| clean(s))
                .collect::<Vec<_>>()
                .join(",")
        ),
    ));
    seq.features.push(feature(
        "misc_feature",
        Location::single(g.upstream_bp as i64),
        "Annotated TSS candidate",
        format!(
            "Genomic {}; annotation-derived, not experimentally established initiation",
            g.tss_1based
        ),
    ));
    for t in &context.transcripts {
        for e in &t.exons {
            let mut f = feature(
                "exon",
                location(&e.span, g),
                &format!("{} E{}", t.transcript_id, e.number_5prime_to_3prime),
                span_note(&e.span),
            );
            f.qualifiers
                .push(("number".into(), Some(e.number_5prime_to_3prime.to_string())));
            f.qualifiers
                .push(("transcript_id".into(), Some(clean(&t.transcript_id))));
            seq.features.push(f);
        }
        for cds in &t.cds {
            seq.features.push(feature(
                "misc_feature",
                location(cds, g),
                &format!("{} CDS segment", t.transcript_id),
                format!(
                    "{}; annotated coding segment; phase not supplied, no translation inferred",
                    span_note(cds)
                ),
            ));
        }
        for c in &t.codons {
            seq.features.push(feature(
                "misc_feature",
                Location::single(c.position_0based as i64),
                &format!("{} translation {:?}", t.transcript_id, c.kind),
                format!("Genomic {}; {}", c.genomic_position_1based, c.basis),
            ));
        }
    }
    for lane in &context.occupancy {
        seq.comments.push(clean(&format!("Signal lane {}/{}: {}; state={}; {} intervals in this window; source_id={}; source_sha256={}; assay={}; factor={}; mark={}; inherited_abs_max={}", lane.group_id, lane.lane_id, lane.label, lane.state.as_str(), lane.intervals.len(), lane.source_id, lane.source_sha256.as_deref().unwrap_or("not supplied"), lane.assay.as_deref().unwrap_or("not supplied"), lane.factor.as_deref().unwrap_or("not supplied"), lane.mark.as_deref().unwrap_or("not supplied"), lane.display_abs_max_score)));
        for interval in &lane.intervals {
            seq.features.push(feature("misc_feature", location(&interval.span, g), &format!("{}: {}", lane.label, interval.interval_id), format!("{}; signal interval, not an individual read; group={}; lane={}; source_id={}; raw_score={}; {}", span_note(&interval.span), lane.group_id, lane.lane_id, lane.source_id, nullable(interval.score), interval.label.as_deref().unwrap_or(""))));
        }
    }
    if let Some(tata) = &context.tata {
        for row in &tata.rows {
            let mut loc = location(&row.span, g);
            if row.genomic_strand != g.strand {
                loc = Location::Complement(Box::new(loc));
            }
            seq.features.push(feature(
                "misc_feature",
                loc,
                &row.evidence.label,
                format!(
                    "{}; evidence_kind={:?}; genomic_strand={}; TATA report={}; report_sha256={}",
                    span_note(&row.span),
                    row.evidence.evidence_kind,
                    row.genomic_strand.as_str(),
                    tata.report_id,
                    tata.report_sha256
                ),
            ));
        }
    }
    for track in &window.tracks {
        for (strand, peaks, maximum) in [
            (
                TssStrand::Plus,
                &track.forward_peaks,
                &track.forward_maximum,
            ),
            (
                TssStrand::Minus,
                &track.reverse_peaks,
                &track.reverse_maximum,
            ),
        ] {
            let mut seen = BTreeSet::new();
            for peak in maximum.iter().chain(peaks) {
                if !seen.insert(peak.local_start_0based) {
                    continue;
                }
                let start = peak.local_start_0based;
                let end = start
                    .checked_add(track.motif_length_bp)
                    .filter(|end| *end <= bases.len())
                    .ok_or_else(|| invalid("GenBank motif span exceeds the TSS window"))?;
                let mut loc = Location::simple_range(start as i64, end as i64);
                if strand == TssStrand::Minus {
                    loc = Location::Complement(Box::new(loc));
                }
                seq.features.push(feature("misc_feature", loc, &format!("{} stored score peak", track.accession), format!("Predicted motif; accession={}; raw_score={}; score_kind={}; motif_local_strand={}; genomic_strand={}; genomic_window_start={}; no new peak calling", track.accession, peak.score, report.panel_resolution.panel.score_kind, strand.as_str(), genomic_strand(g.strand, strand).as_str(), g.genomic_at(start).unwrap())));
            }
        }
    }
    let mut bytes = Vec::new();
    gb_io::writer::write(&mut bytes, &seq)
        .map_err(|e| io_error("write annotated TSS GenBank", e))?;
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Hand-crafted five-base windows using the export-contract fixture; not motif/scorer evidence.
    fn fixture() -> TssProfileReport {
        let mut report = super::super::tests::synthetic_report();
        for w in &mut report.windows {
            let g = &w.record.geometry;
            let span = TssContextSpan {
                genomic_start_1based: g.start_1based - 2,
                genomic_end_1based: g.end_1based,
                start_0based: 0,
                end_0based_exclusive: 5,
                clipped: true,
            };
            w.detail_context = Some(TssDetailContext {
                schema: CONTEXT_SCHEMA.into(),
                window_sequence: Some("ACGTA".into()),
                promoter_id: w.record.promoter_id.clone(),
                geometry: g.clone(),
                window_sequence_sha256: w.record.sequence_sha256.clone(),
                locus_seq_id: "synthetic-locus".into(),
                locus_sequence_sha256: "a".repeat(64),
                locus_report_sha256: "b".repeat(64),
                annotation_release: None,
                bindings: vec![TssInputBinding {
                    role: "tss_detail_locus_report".into(),
                    name: "synthetic.json".into(),
                    sha256: "b".repeat(64),
                }],
                transcripts: vec![TssContextTranscript {
                    transcript_id: w.record.transcripts[0].clone(),
                    label: "Synthetic exon".into(),
                    genomic_strand: g.strand,
                    exons: vec![TssContextExon {
                        number_5prime_to_3prime: 1,
                        span: span.clone(),
                    }],
                    cds: vec![span],
                    codons: vec![],
                }],
                occupancy: vec![],
                tata: None,
                warnings: vec![],
                non_claims: CONTEXT_NON_CLAIMS.into(),
            });
        }
        report
    }

    #[test]
    fn annotated_tss_genbank_roundtrips_strands_clipping_and_raw_peak_scores() {
        let report = fixture();
        let before = serde_json::to_vec(&report).unwrap();
        for window in &report.windows {
            let data = bytes(&report, window).unwrap();
            let seq = gb_io::reader::SeqReader::new(data.as_slice())
                .next()
                .unwrap()
                .unwrap();
            assert_eq!(seq.seq.to_ascii_uppercase(), b"ACGTA");
            let exon = seq.features.iter().find(|f| f.kind == "exon").unwrap();
            assert_eq!(
                exon.location.to_gb_format(),
                if window.record.geometry.strand == TssStrand::Plus {
                    "<1..5"
                } else {
                    "1..>5"
                }
            );
            assert!(!seq.features.iter().any(|f| f.kind == "CDS"));
            let reverse = seq
                .features
                .iter()
                .find(|f| {
                    f.qualifier_values("label")
                        .any(|v| v == "MA0001.1 stored score peak")
                        && matches!(f.location, Location::Complement(_))
                })
                .unwrap();
            assert_eq!(reverse.location.to_gb_format(), "complement(1..3)");
            let note = reverse.qualifier_values("note").next().unwrap();
            assert!(note.contains("raw_score=4"));
            assert!(
                note.contains(if window.record.geometry.strand == TssStrand::Plus {
                    "genomic_strand=-"
                } else {
                    "genomic_strand=+"
                })
            );
            assert_eq!(data, bytes(&report, window).unwrap());
        }
        assert_eq!(before, serde_json::to_vec(&report).unwrap());
    }

    #[test]
    fn annotated_tss_genbank_exports_verify_and_fail_closed_without_bound_bases() {
        let report = fixture();
        let temp = tempfile::tempdir().unwrap();
        let root = fs::canonicalize(temp.path()).unwrap();
        let request = ExportTssProfilesRequest {
            output_dir: root.join("valid").to_str().unwrap().into(),
            context_manifest: None,
            rendering: Default::default(),
            formats: vec![TssExportFormat::Svg, TssExportFormat::Genbank],
        };
        let receipt = export_tss_profiles(&report, &request).unwrap();
        assert_eq!(
            receipt
                .outputs
                .keys()
                .filter(|n| n.ends_with(".gb"))
                .count(),
            2
        );
        read_and_verify_tss_profile_receipt(Path::new(&request.output_dir)).unwrap();
        for (index, bases) in [None, Some("AAAAA".to_string())].into_iter().enumerate() {
            let mut bad = report.clone();
            bad.windows[0]
                .detail_context
                .as_mut()
                .unwrap()
                .window_sequence = bases;
            let request = ExportTssProfilesRequest {
                output_dir: root.join(format!("bad-{index}")).to_str().unwrap().into(),
                ..request.clone()
            };
            assert!(export_tss_profiles(&bad, &request).is_err());
            assert!(!Path::new(&request.output_dir).exists());
        }
        let file = receipt.outputs.keys().find(|n| n.ends_with(".gb")).unwrap();
        fs::write(Path::new(&request.output_dir).join(file), b"tampered").unwrap();
        assert!(read_and_verify_tss_profile_receipt(Path::new(&request.output_dir)).is_err());
    }
}
