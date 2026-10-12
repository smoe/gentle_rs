//! Three-source TATA-box inspection and explicitly approved map annotations.
//! No external software or network access is needed; EPD files are optional.

use super::*;
use gentle_protocol::tata_boxes::*;
use std::io::Read;

fn invalid(message: impl Into<String>) -> EngineError {
    EngineError::new(ErrorCode::InvalidInput, message)
}

fn json_hash(value: &impl serde::Serialize) -> Result<String, EngineError> {
    serde_json::to_vec(value)
        .map(|v| sha256_prefixed_bytes(&v))
        .map_err(|e| invalid(format!("TATA evidence serialization: {e}")))
}

fn exact_location(location: &gb_io::seq::Location) -> bool {
    use gb_io::seq::Location;
    match location {
        Location::Range((start, before), (end, after)) => {
            !before.0 && !after.0 && *start >= 0 && end > start
        }
        Location::Complement(inner) => exact_location(inner),
        Location::Join(parts) => !parts.is_empty() && parts.iter().all(exact_location),
        _ => false,
    }
}

fn text_file(path: &str) -> Result<Option<String>, EngineError> {
    let file = match std::fs::File::open(path) {
        Ok(file) => file,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(e) => return Err(invalid(format!("Cannot read EPD file '{path}': {e}"))),
    };
    let mut bytes = Vec::new();
    file.take(32 * 1024 * 1024 + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| invalid(format!("Cannot read EPD file '{path}': {e}")))?;
    if bytes.len() > 32 * 1024 * 1024 {
        return Err(invalid("EPD input exceeds 32 MiB bound"));
    }
    String::from_utf8(bytes)
        .map(Some)
        .map_err(|e| invalid(format!("EPD input is not UTF-8: {e}")))
}

fn assembly_family(value: &str) -> Option<(&'static str, u64)> {
    let tokens = value
        .split(|c: char| !c.is_ascii_alphanumeric())
        .collect::<Vec<_>>();
    for (family, alias, taxon) in [
        ("grch38", "hg38", 9606),
        ("grch37", "hg19", 9606),
        ("grcm38", "mm10", 10090),
        ("grcm39", "mm39", 10090),
    ] {
        if tokens
            .iter()
            .any(|t| t.eq_ignore_ascii_case(family) || t.eq_ignore_ascii_case(alias))
        {
            return Some((family, taxon));
        }
    }
    None
}

#[derive(Debug)]
struct EpdPromoter {
    chromosome: String,
    tss: usize,
    reverse: bool,
    id: String,
}

fn parse_epd(
    bed: &str,
    motifs: &str,
) -> Result<(Vec<EpdPromoter>, BTreeMap<String, bool>), EngineError> {
    let mut classes = BTreeMap::new();
    let mut lines = motifs.lines().filter(|line| !line.trim().is_empty());
    let header = lines
        .next()
        .ok_or_else(|| invalid("Empty EPD motif table"))?
        .split_whitespace()
        .collect::<Vec<_>>();
    let name_col = header
        .iter()
        .position(|v| *v == "#Name")
        .ok_or_else(|| invalid("EPD motif table requires #Name header"))?;
    let tata_col = header
        .iter()
        .position(|v| *v == "TATA-box")
        .ok_or_else(|| invalid("EPD motif table requires TATA-box header"))?;
    for line in lines {
        let fields = line.split_whitespace().collect::<Vec<_>>();
        if fields.len() != header.len() {
            return Err(invalid("Malformed EPD motif row"));
        }
        let value = match fields[tata_col] {
            "1" => true,
            "0" => false,
            _ => return Err(invalid("EPD TATA-box must be 0 or 1")),
        };
        if classes
            .insert(fields[name_col].to_string(), value)
            .is_some()
        {
            return Err(invalid("Duplicate EPD motif promoter ID"));
        }
    }
    let mut promoters = Vec::new();
    let mut ids = BTreeSet::new();
    for (line_index, line) in bed.lines().enumerate() {
        let line = line.trim();
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track ")
            || line.starts_with("browser ")
        {
            continue;
        }
        let fields = line.split_whitespace().collect::<Vec<_>>();
        let error = || {
            invalid(format!(
                "Invalid EPD BED8 line {}: expected official 60-bp promoter / 11-bp thick TSS segment",
                line_index + 1
            ))
        };
        if fields.len() < 8 {
            return Err(error());
        }
        let number = |i: usize| fields[i].parse::<usize>().map_err(|_| error());
        let (start, end, thick_start, thick_end) = (number(1)?, number(2)?, number(6)?, number(7)?);
        let reverse = match fields[5] {
            "+" => false,
            "-" => true,
            _ => return Err(error()),
        };
        if end.checked_sub(start) != Some(60)
            || thick_end.checked_sub(thick_start) != Some(11)
            || thick_start < start
            || thick_end > end
            || (reverse && thick_start != start)
            || (!reverse && thick_end != end)
        {
            return Err(error());
        }
        if !ids.insert(fields[3].to_string()) {
            return Err(invalid("Duplicate EPD BED promoter ID"));
        }
        promoters.push(EpdPromoter {
            chromosome: fields[0].to_string(),
            tss: if reverse { thick_end - 1 } else { thick_start },
            reverse,
            id: fields[3].to_string(),
        });
    }
    if promoters.is_empty() {
        return Err(invalid("EPD BED has no promoter records"));
    }
    Ok((promoters, classes))
}

fn associated_tss(
    start: usize,
    end: usize,
    reverse: bool,
    tss: &[TataBoxTss],
    request: &TataBoxScreenRequest,
) -> Vec<TataBoxTssAssociation> {
    let first = if reverse { end - 1 } else { start } as i64;
    tss.iter()
        .filter(|t| t.reverse == reverse)
        .filter_map(|t| {
            let distance = if reverse {
                t.position_0based as i64 - first
            } else {
                first - t.position_0based as i64
            };
            (distance >= i64::from(request.minimum_tss_offset)
                && distance <= i64::from(request.maximum_tss_offset))
            .then(|| TataBoxTssAssociation {
                tss_id: t.id.clone(),
                source: t.source.clone(),
                signed_distance_bp: distance,
            })
        })
        .collect()
}

fn search_ranges(
    request: &TataBoxScreenRequest,
    end: usize,
    width: usize,
    reverse: bool,
    tss: &[TataBoxTss],
) -> Vec<(usize, usize)> {
    if width > end - request.start_0based {
        return vec![];
    }
    let first = request.start_0based as i64;
    let last_exclusive = (end - width + 1) as i64;
    if request.scan_without_tss {
        return vec![(first as usize, last_exclusive as usize)];
    }
    let mut ranges = tss
        .iter()
        .filter(|t| t.reverse == reverse)
        .filter_map(|t| {
            let position = t.position_0based as i64;
            let (start, stop) = if reverse {
                (
                    position - i64::from(request.maximum_tss_offset) - width as i64 + 1,
                    position - i64::from(request.minimum_tss_offset) - width as i64 + 2,
                )
            } else {
                (
                    position + i64::from(request.minimum_tss_offset),
                    position + i64::from(request.maximum_tss_offset) + 1,
                )
            };
            let (start, stop) = (start.max(first), stop.min(last_exclusive));
            (start < stop).then_some((start as usize, stop as usize))
        })
        .collect::<Vec<_>>();
    ranges.sort_unstable();
    let mut merged: Vec<(usize, usize)> = vec![];
    for (start, stop) in ranges {
        if let Some(last) = merged.last_mut()
            && start <= last.1
        {
            last.1 = last.1.max(stop);
        } else {
            merged.push((start, stop));
        }
    }
    merged
}

impl GentleEngine {
    /// Deterministic read-only inspection; no missing input triggers a download.
    pub fn screen_tata_boxes(
        &self,
        request: &TataBoxScreenRequest,
    ) -> Result<TataBoxScreenReport, EngineError> {
        let dna = self
            .state
            .sequences
            .get(&request.seq_id)
            .ok_or_else(|| invalid(format!("Unknown sequence '{}'", request.seq_id)))?;
        let end = request.end_0based_exclusive.unwrap_or(dna.len());
        if request.start_0based >= end
            || end > dna.len()
            || end - request.start_0based > 2_000_000
            || dna.len() > i64::MAX as usize
            || request.max_rows == 0
            || request.max_rows > 50_000
            || !request.minimum_llr_bits.is_finite()
            || request.minimum_tss_offset > request.maximum_tss_offset
            || request.minimum_tss_offset < -100_000
            || request.maximum_tss_offset > 100_000
        {
            return Err(invalid(
                "TATA screen needs a valid range (at most 2 Mb), finite score, ordered offsets within +/-100000, and 1..50000 rows",
            ));
        }
        let mut effective_request = request.clone();
        effective_request.end_0based_exclusive = Some(end);
        let anchor = self.sequence_genome_anchor_summary(&request.seq_id).ok();
        let mut report = TataBoxScreenReport {
            schema: TATA_BOX_SCREEN_SCHEMA.into(), report_id: String::new(), content_sha256: String::new(),
            request: effective_request, sequence_sha256: sha256_prefixed_bytes(dna.forward_bytes()),
            annotation_sha256: json_hash(dna.features())?, anchor_sha256: json_hash(&anchor)?,
            genome_id: anchor.as_ref().map(|a| a.genome_id.clone()), chromosome: anchor.as_ref().map(|a| a.chromosome.clone()),
            genomic_start_1based: anchor.as_ref().map(|a| a.start_1based), genomic_end_1based: anchor.as_ref().map(|a| a.end_1based),
            genomic_reverse: anchor.as_ref().map(|a| a.strand == Some('-')),
            motif_id: None, matrix_sha256: None,
            score_policy: "GENtle LLR bits; shared smoothing; uniform A/C/G/T background; threshold is not a p-value".into(),
            epd_status: "not_requested".into(), epd_bed_sha256: None, epd_motifs_sha256: None,
            tss: vec![], rows: vec![], scored_windows: 0, ambiguous_windows: 0, warnings: vec![], non_claim: TATA_BOX_NON_CLAIM.into(),
        };
        for t in &request.additional_tss {
            if t.id.trim().is_empty()
                || t.source.trim().is_empty()
                || t.position_0based >= dna.len()
            {
                return Err(invalid(
                    "Explicit TSS requires id, source and an in-sequence position",
                ));
            }
            report.tss.push(t.clone());
        }
        for t in self.derive_promoter_window_records(
            dna,
            None,
            None,
            0,
            0,
            PromoterWindowCollapseMode::Transcript,
        ) {
            if t.transcript_feature_id
                .and_then(|id| dna.features().get(id))
                .is_some_and(|f| exact_location(&f.location))
            {
                report.tss.push(TataBoxTss {
                    id: t.transcript_id,
                    position_0based: t.tss_local_0based,
                    reverse: t.strand == "-",
                    source: "annotated_transcript_5prime_boundary_not_experimentally_confirmed_TSS"
                        .into(),
                });
            } else {
                report.warnings.push(format!(
                    "Skipped imprecise transcript boundary for {}",
                    t.transcript_id
                ));
            }
        }
        if let Some(source) = &request.epd {
            self.tata_epd_rows(source, dna, anchor.as_ref(), &mut report)?;
        }
        report.tss.sort();
        report.tss.dedup();
        if request.include_annotations {
            for (feature_id, feature) in dna.features().iter().enumerate() {
                if feature
                    .qualifier_values("gentle_generated")
                    .any(|v| v == "tata_box_screen")
                {
                    continue;
                }
                let is_tata = feature.kind.eq_ignore_ascii_case("TATA_signal")
                    || (feature.kind.eq_ignore_ascii_case("regulatory")
                        && feature
                            .qualifier_values("regulatory_class")
                            .any(|v| v.eq_ignore_ascii_case("TATA_box")));
                if !is_tata {
                    continue;
                }
                let mut ranges = vec![];
                collect_location_ranges_usize(&feature.location, &mut ranges);
                if !exact_location(&feature.location) || ranges.len() != 1 {
                    report.warnings.push(format!("TATA annotation feature {feature_id} has imprecise/compound geometry; not projected as an exact site"));
                    continue;
                }
                let (start, stop) = ranges[0];
                if start < request.start_0based || stop > end {
                    continue;
                }
                let reverse = feature_is_reverse(feature);
                report.rows.push(TataBoxEvidenceRow {
                    row_id: format!("annotation_{feature_id}"),
                    evidence_kind: TataBoxEvidenceKind::SourceAnnotation,
                    label: Self::feature_display_label(feature, feature_id),
                    start_0based: start,
                    end_0based_exclusive: stop,
                    reverse,
                    geometry_kind: "exact_site".into(),
                    source_feature_id: Some(feature_id),
                    source_qualifiers: feature
                        .qualifiers
                        .iter()
                        .map(|(k, v)| (k.to_string(), v.clone()))
                        .collect(),
                    tata_positive: None,
                    llr_bits: None,
                    sequence_5prime_to_3prime: None,
                    tss_associations: associated_tss(start, stop, reverse, &report.tss, request),
                });
            }
        }
        if request.predict {
            let motif = crate::tf_motifs::resolve_motif_definition(&request.motif_id)
                .ok_or_else(|| invalid(format!("TBP matrix '{}' unavailable", request.motif_id)))?;
            if motif.id != request.motif_id
                || !motif.has_full_pfm
                || !motif
                    .name
                    .as_deref()
                    .is_some_and(|v| v.eq_ignore_ascii_case("TBP"))
                || motif.matrix_counts.is_empty()
                || motif.matrix_counts.iter().any(|column| {
                    column.iter().any(|v| !v.is_finite() || *v < 0.0)
                        || column.iter().sum::<f64>() <= 0.0
                })
            {
                return Err(invalid(
                    "TATA prediction requires an exact versioned full-PFM TBP matrix ID",
                ));
            }
            let (llr, _) = Self::prepare_scoring_matrices(&motif.matrix_counts);
            report.motif_id = Some(motif.id.clone());
            report.matrix_sha256 = Some(json_hash(&motif.matrix_counts)?);
            let width = llr.len();
            if width <= end - request.start_0based {
                let sequence = dna.forward_bytes();
                for reverse in [false, true] {
                    for start in search_ranges(request, end, width, reverse, &report.tss)
                        .into_iter()
                        .flat_map(|(start, stop)| start..stop)
                    {
                        let forward = &sequence[start..start + width];
                        let reversed;
                        let window = if reverse {
                            reversed = forward
                                .iter()
                                .rev()
                                .map(|b| match b.to_ascii_uppercase() {
                                    b'A' => b'T',
                                    b'C' => b'G',
                                    b'G' => b'C',
                                    b'T' => b'A',
                                    _ => b'N',
                                })
                                .collect::<Vec<_>>();
                            reversed.as_slice()
                        } else {
                            forward
                        };
                        let Some(score) = Self::score_matrix_window(window, &llr) else {
                            report.ambiguous_windows += 1;
                            continue;
                        };
                        report.scored_windows += 1;
                        if score < request.minimum_llr_bits {
                            continue;
                        }
                        let associations =
                            associated_tss(start, start + width, reverse, &report.tss, request);
                        report.rows.push(TataBoxEvidenceRow {
                            row_id: format!(
                                "tbp_{}_{}_{}",
                                motif.id,
                                start,
                                if reverse { "minus" } else { "plus" }
                            ),
                            evidence_kind: TataBoxEvidenceKind::MotifPrediction,
                            label: format!("TBP {} candidate", motif.id),
                            start_0based: start,
                            end_0based_exclusive: start + width,
                            reverse,
                            geometry_kind: "exact_site".into(),
                            source_feature_id: None,
                            source_qualifiers: vec![],
                            tata_positive: None,
                            llr_bits: Some(score),
                            sequence_5prime_to_3prime: Some(
                                String::from_utf8_lossy(window).to_ascii_uppercase(),
                            ),
                            tss_associations: associations,
                        });
                        if report.rows.len() > request.max_rows {
                            return Err(invalid(
                                "TATA result exceeds max_rows; narrow the range or increase the score threshold",
                            ));
                        }
                    }
                }
            }
            if !request.scan_without_tss && report.tss.is_empty() {
                report.warnings.push("No exact TSS context available; motif prediction not evaluated. Supply a TSS or explicitly enable scan_without_tss.".into());
            }
        }
        if report.rows.len() > request.max_rows {
            return Err(invalid("TATA evidence exceeds max_rows; narrow the range"));
        }
        report.rows.sort_by(|a, b| {
            (
                a.start_0based,
                a.end_0based_exclusive,
                a.reverse,
                a.evidence_kind,
                &a.row_id,
            )
                .cmp(&(
                    b.start_0based,
                    b.end_0based_exclusive,
                    b.reverse,
                    b.evidence_kind,
                    &b.row_id,
                ))
        });
        report.warnings.sort();
        report.warnings.dedup();
        report.report_id = short_sha256_id("tata_screen", &json_hash(&report)?);
        report.content_sha256 = json_hash(&report)?;
        Ok(report)
    }

    fn tata_epd_rows(
        &self,
        source: &TataBoxEpdSource,
        dna: &DNAsequence,
        anchor: Option<&SequenceGenomeAnchorSummary>,
        report: &mut TataBoxScreenReport,
    ) -> Result<(), EngineError> {
        if source.release.trim().is_empty()
            || source.source_url.trim().is_empty()
            || source.bed_path.trim().is_empty()
            || source.motifs_path.trim().is_empty()
        {
            return Err(invalid(
                "EPD requires release, source_url, bed_path and motifs_path",
            ));
        }
        let Some(anchor) = anchor else {
            if source.required {
                return Err(invalid(
                    "Required EPD source needs a genome-anchored sequence",
                ));
            }
            report.epd_status = "unavailable_no_genome_anchor".into();
            return Ok(());
        };
        let family = assembly_family(&source.assembly).ok_or_else(|| {
            invalid("EPD assembly must be GRCh37/hg19, GRCh38/hg38, GRCm38/mm10 or GRCm39/mm39")
        })?;
        if assembly_family(&anchor.genome_id) != Some(family) || family.1 != source.taxon_id {
            return Err(invalid(format!(
                "EPD assembly/species mismatch: {} taxon {} versus sequence genome {}",
                source.assembly, source.taxon_id, anchor.genome_id
            )));
        }
        if anchor.strand.is_none()
            || anchor.start_1based == 0
            || anchor
                .end_1based
                .checked_sub(anchor.start_1based)
                .and_then(|v| v.checked_add(1))
                != Some(dna.len())
        {
            return Err(invalid(
                "EPD projection requires an exact full-length strand-aware genome anchor",
            ));
        }
        let bed = text_file(&source.bed_path)?;
        let motifs = text_file(&source.motifs_path)?;
        let (Some(bed), Some(motifs)) = (bed, motifs) else {
            if source.required {
                return Err(invalid("Required EPD BED or motif table is missing"));
            }
            report.epd_status = "unavailable_missing_files".into();
            report
                .warnings
                .push("EPD files missing; no EPD classification was inferred".into());
            return Ok(());
        };
        let bed_hash = sha256_prefixed_str(&bed);
        let motif_hash = sha256_prefixed_str(&motifs);
        if source
            .expected_bed_sha256
            .as_ref()
            .is_some_and(|v| *v != bed_hash)
            || source
                .expected_motifs_sha256
                .as_ref()
                .is_some_and(|v| *v != motif_hash)
        {
            return Err(invalid("EPD input content hash mismatch"));
        }
        report.epd_bed_sha256 = Some(bed_hash);
        report.epd_motifs_sha256 = Some(motif_hash);
        let (promoters, classes) = parse_epd(&bed, &motifs)?;
        report.epd_status = "available".into();
        for promoter in promoters {
            if !Self::genome_chromosome_matches(&promoter.chromosome, &anchor.chromosome)
                || promoter.tss < anchor.start_1based - 1
                || promoter.tss >= anchor.end_1based
            {
                continue;
            }
            let local = if anchor.strand == Some('-') {
                anchor.end_1based - promoter.tss - 1
            } else {
                promoter.tss - (anchor.start_1based - 1)
            };
            let reverse = promoter.reverse != (anchor.strand == Some('-'));
            report.tss.push(TataBoxTss {
                id: promoter.id.clone(),
                position_0based: local,
                reverse,
                source: format!("EPDnew:{}", source.release),
            });
            if local < report.request.start_0based
                || local >= report.request.end_0based_exclusive.unwrap_or(dna.len())
            {
                continue;
            }
            let positive = classes.get(&promoter.id).copied();
            if positive.is_none() {
                report.warnings.push(format!(
                    "EPD promoter {} has no TATA classification row",
                    promoter.id
                ));
            }
            report.rows.push(TataBoxEvidenceRow {
                row_id: format!("epd_{}", promoter.id),
                evidence_kind: TataBoxEvidenceKind::EpdClassification,
                label: promoter.id,
                start_0based: local,
                end_0based_exclusive: local + 1,
                reverse,
                geometry_kind: "tss_marker_not_tata_site".into(),
                source_feature_id: None,
                source_qualifiers: vec![],
                tata_positive: positive,
                llr_bits: None,
                sequence_5prime_to_3prime: None,
                tss_associations: vec![],
            });
        }
        Ok(())
    }

    /// Add only reviewed rows, retaining source identity. Original annotations
    /// are already on the map and are not duplicated; EPD rows become TSS markers.
    pub fn materialize_tata_boxes(
        &mut self,
        request: &TataBoxMaterializeRequest,
    ) -> Result<TataBoxScreenReport, EngineError> {
        let report = self.screen_tata_boxes(&request.screen)?;
        if report.content_sha256 != request.expected_report_sha256 {
            return Err(invalid(
                "TATA evidence changed since review; inspect again before adding features",
            ));
        }
        let selected = request.row_ids.iter().collect::<BTreeSet<_>>();
        if selected.is_empty()
            || selected.len() != request.row_ids.len()
            || selected
                .iter()
                .any(|id| !report.rows.iter().any(|r| &r.row_id == *id))
        {
            return Err(invalid("Select unique existing TATA evidence row IDs"));
        }
        let dna = self
            .state
            .sequences
            .get_mut(&request.screen.seq_id)
            .ok_or_else(|| invalid("TATA sequence no longer available"))?;
        let original_feature_count = dna.features().len();
        for row in report.rows.iter().filter(|r| selected.contains(&r.row_id)) {
            if row.evidence_kind == TataBoxEvidenceKind::SourceAnnotation {
                continue;
            }
            let key = json_hash(&(
                &report.sequence_sha256,
                &report.anchor_sha256,
                row,
                &report.matrix_sha256,
                &report.epd_bed_sha256,
                &report.epd_motifs_sha256,
            ))?;
            if dna
                .features()
                .iter()
                .any(|f| f.qualifier_values("tata_evidence_id").any(|v| v == key))
            {
                continue;
            }
            let mut location = gb_io::seq::Location::simple_range(
                row.start_0based as i64,
                row.end_0based_exclusive as i64,
            );
            if row.reverse {
                location = gb_io::seq::Location::Complement(Box::new(location));
            }
            let classification = match row.tata_positive {
                Some(true) => "positive",
                Some(false) => "negative",
                None => "unavailable",
            };
            let label = if row.evidence_kind == TataBoxEvidenceKind::EpdClassification {
                format!("EPD {}: TATA {classification} (TSS)", row.label)
            } else {
                format!("TATA candidate {}", row.label)
            };
            dna.features_mut().push(gb_io::seq::Feature {
                kind: "regulatory".into(),
                location,
                qualifiers: vec![
                    ("label".into(), Some(label)),
                    (
                        "regulatory_class".into(),
                        Some(
                            if row.evidence_kind == TataBoxEvidenceKind::EpdClassification {
                                "promoter"
                            } else {
                                "TATA_box"
                            }
                            .into(),
                        ),
                    ),
                    ("gentle_generated".into(), Some("tata_box_screen".into())),
                    (
                        "gentle_feature_source".into(),
                        Some(
                            if row.evidence_kind == TataBoxEvidenceKind::EpdClassification {
                                "EPDnew_classification"
                            } else {
                                "GENtle_TBP_prediction"
                            }
                            .into(),
                        ),
                    ),
                    ("tata_evidence_id".into(), Some(key)),
                    (
                        "tata_report_sha256".into(),
                        Some(report.content_sha256.clone()),
                    ),
                    (
                        "tata_evidence_json".into(),
                        Some(serde_json::to_string(row).map_err(|e| invalid(e.to_string()))?),
                    ),
                    (
                        "tata_request_json".into(),
                        Some(
                            serde_json::to_string(&report.request)
                                .map_err(|e| invalid(e.to_string()))?,
                        ),
                    ),
                    ("tata_matrix_sha256".into(), report.matrix_sha256.clone()),
                    ("tata_epd_bed_sha256".into(), report.epd_bed_sha256.clone()),
                    (
                        "tata_epd_motifs_sha256".into(),
                        report.epd_motifs_sha256.clone(),
                    ),
                    ("note".into(), Some(TATA_BOX_NON_CLAIM.into())),
                ],
            });
        }
        if dna.features().len() != original_feature_count {
            Self::prepare_sequence(dna);
        }
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gb_io::seq::{Feature, Location};
    use tempfile::tempdir;

    // Hand-crafted, non-biological test sequences and EPD-format rows. No
    // experimental promoter or vendor data are copied. Recreated by these tests.
    fn engine(reverse: bool) -> GentleEngine {
        let mut bases = vec![b'C'; 220];
        bases[72..79].copy_from_slice(b"TATAAAA");
        let bases = if reverse {
            bases
                .iter()
                .rev()
                .map(|b| match b {
                    b'C' => b'G',
                    b'T' => b'A',
                    _ => b'T',
                })
                .collect::<Vec<_>>()
        } else {
            bases
        };
        let mut dna = DNAsequence::from_sequence(&String::from_utf8(bases).unwrap()).unwrap();
        let location = if reverse {
            Location::Complement(Box::new(Location::simple_range(20, 120)))
        } else {
            Location::simple_range(100, 200)
        };
        dna.features_mut().push(Feature {
            kind: "mRNA".into(),
            location,
            qualifiers: vec![("transcript_id".into(), Some("synthetic_transcript".into()))],
        });
        let mut engine = GentleEngine::new();
        engine.state_mut().sequences.insert("toy".into(), dna);
        engine.state_mut().metadata.insert(PROVENANCE_METADATA_KEY.into(), serde_json::json!({GENOME_EXTRACTIONS_METADATA_KEY: [{
            "seq_id":"toy", "genome_id":"GRCh38.p14", "chromosome":"1", "start_1based":1001, "end_1based":1220, "anchor_strand":if reverse {"-"} else {"+"}
        }]}));
        GentleEngine::from_state(engine.snapshot().clone())
    }

    fn request() -> TataBoxScreenRequest {
        TataBoxScreenRequest {
            seq_id: "toy".into(),
            ..Default::default()
        }
    }

    fn epd_source(directory: &std::path::Path) -> TataBoxEpdSource {
        let bed_path = directory.join("toy.bed");
        let motifs_path = directory.join("motifs.txt");
        std::fs::write(&bed_path, "chr1 1051 1111 SYNTH_1 900 + 1100 1111\nchr1 1150 1210 SYNTH_2 900 - 1150 1161\nchr1 1061 1121 SYNTH_3 900 + 1110 1121\n").unwrap();
        std::fs::write(
            &motifs_path,
            "#Name\tTATA-box\tInr\nSYNTH_1\t1\t0\nSYNTH_2\t0\t1\n",
        )
        .unwrap();
        TataBoxEpdSource {
            bed_path: bed_path.to_string_lossy().into(),
            motifs_path: motifs_path.to_string_lossy().into(),
            assembly: "GRCh38".into(),
            taxon_id: 9606,
            release: "synthetic-001".into(),
            source_url: "synthetic:epd-format".into(),
            required: true,
            expected_bed_sha256: None,
            expected_motifs_sha256: None,
        }
    }

    #[test]
    fn tata_tss_scan_is_strand_aware_and_deterministic_read_only() {
        let _lock = crate::tf_motifs::test_registry_lock().lock().unwrap();
        for reverse in [false, true] {
            let mut engine = engine(reverse);
            let before = serde_json::to_value(engine.state()).unwrap();
            let a = engine.screen_tata_boxes(&request()).unwrap();
            let b = engine.screen_tata_boxes(&request()).unwrap();
            assert_eq!(a, b);
            let row = a
                .rows
                .iter()
                .find(|r| r.sequence_5prime_to_3prime.as_deref() == Some("TATAAAA"))
                .unwrap();
            assert_eq!(row.reverse, reverse);
            assert_eq!(
                (row.start_0based, row.end_0based_exclusive),
                if reverse { (141, 148) } else { (72, 79) }
            );
            assert_eq!(row.tss_associations[0].signed_distance_bp, -28);
            assert_eq!(a.scored_windows, 26);
            assert_eq!(a.motif_id.as_deref(), Some("MA0108.3"));
            assert!(a.matrix_sha256.is_some());
            let mut hashable = a.clone();
            hashable.content_sha256.clear();
            assert_eq!(a.content_sha256, json_hash(&hashable).unwrap());
            let result = engine
                .apply(Operation::ScreenTataBoxes {
                    request: request(),
                    path: None,
                })
                .unwrap();
            assert!(result.tata_box_screen.is_some());
            assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
            assert_eq!(engine.undo_available(), 0);
        }
    }

    #[test]
    fn tata_source_annotations_are_separate_and_keep_provenance() {
        let _lock = crate::tf_motifs::test_registry_lock().lock().unwrap();
        let mut engine = engine(false);
        let dna = engine.state_mut().sequences.get_mut("toy").unwrap();
        dna.features_mut().push(Feature {
            kind: "regulatory".into(),
            location: Location::simple_range(72, 79),
            qualifiers: vec![
                ("regulatory_class".into(), Some("TATA_box".into())),
                ("inference".into(), Some("synthetic:test".into())),
            ],
        });
        dna.features_mut().push(Feature {
            kind: "TATA_signal".into(),
            location: Location::simple_range(40, 47),
            qualifiers: vec![],
        });
        dna.features_mut().push(Feature {
            kind: "misc_feature".into(),
            location: Location::simple_range(20, 27),
            qualifiers: vec![("note".into(), Some("not a curated TATA annotation".into()))],
        });
        let report = engine.screen_tata_boxes(&request()).unwrap();
        assert_eq!(
            report
                .rows
                .iter()
                .filter(|r| r.evidence_kind == TataBoxEvidenceKind::SourceAnnotation)
                .count(),
            2
        );
        let row = report
            .rows
            .iter()
            .find(|r| r.source_feature_id == Some(1))
            .unwrap();
        assert!(
            row.source_qualifiers
                .iter()
                .any(|(k, v)| k == "inference" && v.as_deref() == Some("synthetic:test"))
        );
        assert_eq!(row.llr_bits, None);
    }

    #[test]
    fn tata_epd_bed8_tss_classification_never_invents_site_coordinates() {
        let _lock = crate::tf_motifs::test_registry_lock().lock().unwrap();
        let dir = tempdir().unwrap();
        let source = epd_source(dir.path());
        for reverse in [false, true] {
            let engine = engine(reverse);
            let mut request = request();
            request.epd = Some(source.clone());
            let report = engine.screen_tata_boxes(&request).unwrap();
            let epd = report
                .rows
                .iter()
                .filter(|r| r.evidence_kind == TataBoxEvidenceKind::EpdClassification)
                .collect::<Vec<_>>();
            assert_eq!(epd.len(), 3);
            let positive = epd.iter().find(|r| r.label == "SYNTH_1").unwrap();
            assert_eq!(positive.start_0based, if reverse { 119 } else { 100 });
            assert_eq!(positive.end_0based_exclusive - positive.start_0based, 1);
            assert_eq!(positive.geometry_kind, "tss_marker_not_tata_site");
            assert_eq!(positive.reverse, reverse);
            assert_eq!(positive.tata_positive, Some(true));
            let negative = epd.iter().find(|r| r.label == "SYNTH_2").unwrap();
            assert_eq!(negative.start_0based, if reverse { 59 } else { 160 });
            assert_eq!(negative.reverse, !reverse);
            assert_eq!(negative.tata_positive, Some(false));
            assert_eq!(
                epd.iter()
                    .find(|r| r.label == "SYNTH_3")
                    .unwrap()
                    .tata_positive,
                None
            );
            assert!(report.epd_bed_sha256.is_some() && report.epd_motifs_sha256.is_some());
        }
    }

    #[test]
    fn tata_epd_missing_mismatched_and_stale_sources_fail_honestly() {
        let _lock = crate::tf_motifs::test_registry_lock().lock().unwrap();
        let dir = tempdir().unwrap();
        let source = epd_source(dir.path());
        let engine = engine(false);
        let mut req = request();
        req.epd = Some(source.clone());
        for wrong in ["assembly", "taxon", "hash"] {
            let mut bad = source.clone();
            match wrong {
                "assembly" => bad.assembly = "hg19".into(),
                "taxon" => bad.taxon_id = 10090,
                _ => bad.expected_motifs_sha256 = Some("sha256:wrong".into()),
            }
            req.epd = Some(bad);
            assert!(engine.screen_tata_boxes(&req).is_err());
        }
        std::fs::remove_file(&source.bed_path).unwrap();
        req.epd = Some(source.clone());
        assert!(engine.screen_tata_boxes(&req).is_err());
        req.epd.as_mut().unwrap().required = false;
        let report = engine.screen_tata_boxes(&req).unwrap();
        assert_eq!(report.epd_status, "unavailable_missing_files");
        assert!(
            report
                .rows
                .iter()
                .any(|r| r.evidence_kind == TataBoxEvidenceKind::MotifPrediction)
        );
        assert!(
            !report
                .rows
                .iter()
                .any(|r| r.evidence_kind == TataBoxEvidenceKind::EpdClassification)
        );
    }

    #[test]
    fn tata_epd_parser_rejects_generic_bed_and_duplicate_ids() {
        assert!(parse_epd("chr1 10 20 X 0 +\n", "#Name TATA-box\nX 1\n").is_err());
        let bed = "chr1 51 111 X 900 + 100 111\n";
        assert!(parse_epd(bed, "#Name TATA-box\nX 1\nX 0\n").is_err());
        assert!(parse_epd(bed, "#Name TATA-box\nX NA\n").is_err());
    }

    #[test]
    fn tata_missing_tss_and_ambiguous_sequence_are_not_negative_evidence() {
        let _lock = crate::tf_motifs::test_registry_lock().lock().unwrap();
        let mut engine = engine(false);
        engine
            .state_mut()
            .sequences
            .get_mut("toy")
            .unwrap()
            .features_mut()
            .clear();
        let report = engine.screen_tata_boxes(&request()).unwrap();
        assert_eq!(report.scored_windows, 0);
        assert!(!report.warnings.is_empty());
        let mut req = request();
        req.scan_without_tss = true;
        assert!(engine.screen_tata_boxes(&req).unwrap().scored_windows > 0);
        engine.state_mut().sequences.insert(
            "toy".into(),
            DNAsequence::from_sequence(&"N".repeat(100)).unwrap(),
        );
        let report = engine.screen_tata_boxes(&req).unwrap();
        assert_eq!(report.scored_windows, 0);
        assert_eq!(report.ambiguous_windows, 188);
        req.max_rows = 1;
        req.minimum_llr_bits = -1000.0;
        engine.state_mut().sequences.insert(
            "toy".into(),
            DNAsequence::from_sequence(&"A".repeat(100)).unwrap(),
        );
        assert!(
            engine
                .screen_tata_boxes(&req)
                .unwrap_err()
                .message
                .contains("max_rows")
        );
    }

    #[test]
    fn tata_materialization_is_digest_gated_explicit_and_undoable() {
        let _lock = crate::tf_motifs::test_registry_lock().lock().unwrap();
        let dir = tempdir().unwrap();
        let mut engine = engine(false);
        let mut req = request();
        req.epd = Some(epd_source(dir.path()));
        let report = engine.screen_tata_boxes(&req).unwrap();
        let before = serde_json::to_value(engine.state()).unwrap();
        let mut apply = TataBoxMaterializeRequest {
            screen: req,
            expected_report_sha256: "wrong".into(),
            row_ids: report.rows.iter().map(|r| r.row_id.clone()).collect(),
        };
        assert!(
            engine
                .apply(Operation::MaterializeTataBoxFeatures {
                    request: apply.clone()
                })
                .is_err()
        );
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
        apply.expected_report_sha256 = report.content_sha256;
        engine
            .apply(Operation::MaterializeTataBoxFeatures {
                request: apply.clone(),
            })
            .unwrap();
        let features = engine.state().sequences["toy"].features();
        assert!(features.iter().any(|f| {
            f.qualifier_values("regulatory_class")
                .any(|v| v == "TATA_box")
        }));
        assert!(features.iter().any(|f| {
            f.qualifier_values("gentle_feature_source")
                .any(|v| v == "EPDnew_classification")
                && f.qualifier_values("regulatory_class")
                    .any(|v| v == "promoter")
        }));
        assert!(
            engine
                .apply(Operation::MaterializeTataBoxFeatures { request: apply })
                .is_err()
        );
        engine.undo_last_operation().unwrap();
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    }

    #[test]
    fn tata_reinspection_does_not_duplicate_map_features() {
        let _lock = crate::tf_motifs::test_registry_lock().lock().unwrap();
        let mut engine = engine(false);
        for _ in 0..2 {
            let report = engine.screen_tata_boxes(&request()).unwrap();
            engine
                .apply(Operation::MaterializeTataBoxFeatures {
                    request: TataBoxMaterializeRequest {
                        screen: request(),
                        expected_report_sha256: report.content_sha256,
                        row_ids: report.rows.iter().map(|r| r.row_id.clone()).collect(),
                    },
                })
                .unwrap();
        }
        let count = engine.state().sequences["toy"]
            .features()
            .iter()
            .filter(|f| {
                f.qualifier_values("gentle_generated")
                    .any(|v| v == "tata_box_screen")
            })
            .count();
        assert_eq!(
            count,
            engine.screen_tata_boxes(&request()).unwrap().rows.len()
        );
    }

    #[test]
    fn tata_epd_published_files_opt_in() {
        let Ok(bed_path) = std::env::var("GENTLE_TEST_EPD_BED") else {
            return;
        };
        let motifs_path = std::env::var("GENTLE_TEST_EPD_MOTIFS")
            .expect("GENTLE_TEST_EPD_MOTIFS accompanies EPD BED");
        let (rows, classes) = parse_epd(
            &std::fs::read_to_string(bed_path).unwrap(),
            &std::fs::read_to_string(motifs_path).unwrap(),
        )
        .unwrap();
        assert!(rows.len() > 1000);
        assert!(classes.values().any(|v| *v) && classes.values().any(|v| !v));
        assert!(rows.iter().all(|r| classes.contains_key(&r.id)));
    }
}
