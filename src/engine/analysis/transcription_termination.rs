//! Shared, provenance-bound termination recognition and reviewed map annotation.
//!
//! Searches inspect source DNA, never infer mature transcript ends, and never
//! invoke a folding tool. Explicit feature context restricts strand; intervals
//! may deliberately extend downstream of that feature. Wrapped spans are refused.

use super::*;
use base64::{Engine as _, engine::general_purpose::STANDARD};
use gentle_protocol::transcription_termination::*;
use std::io::Read;

const MAX_CATALOG_BYTES: usize = 4 * 1024 * 1024;
fn invalid(message: impl Into<String>) -> EngineError {
    EngineError::new(ErrorCode::InvalidInput, message)
}
fn hash(value: &impl serde::Serialize) -> Result<String, EngineError> {
    serde_json::to_vec(value)
        .map(|v| sha256_prefixed_bytes(&v))
        .map_err(|e| invalid(format!("Termination report serialization: {e}")))
}
/// GenBank may fold quoted qualifiers at arbitrary positions. Base64 protects
/// JSON string content; only transport folding whitespace is discarded on decode.
fn encode_provenance(value: &impl serde::Serialize) -> Result<String, EngineError> {
    serde_json::to_vec(value)
        .map(|bytes| STANDARD.encode(bytes))
        .map_err(|e| invalid(e.to_string()))
}
fn decode_provenance<T: serde::de::DeserializeOwned>(text: &str) -> Result<T, EngineError> {
    let compact = text
        .bytes()
        .filter(|b| !b.is_ascii_whitespace())
        .collect::<Vec<_>>();
    let bytes = STANDARD
        .decode(compact)
        .map_err(|e| invalid(format!("Invalid stored termination base64: {e}")))?;
    serde_json::from_slice(&bytes)
        .map_err(|e| invalid(format!("Invalid stored termination JSON: {e}")))
}
fn compact_digest(text: &str) -> String {
    text.chars().filter(|c| !c.is_ascii_whitespace()).collect()
}

fn exact_location(location: &gb_io::seq::Location) -> bool {
    use gb_io::seq::Location;
    match location {
        Location::Range((s, before), (e, after)) => !before.0 && !after.0 && *s >= 0 && e > s,
        Location::Complement(inner) => exact_location(inner),
        Location::Join(parts) => !parts.is_empty() && parts.iter().all(exact_location),
        _ => false,
    }
}
fn annotation_geometry(feature: &gb_io::seq::Feature) -> Option<TerminationGeometryKind> {
    if feature.kind.eq_ignore_ascii_case("terminator")
        || feature
            .qualifier_values("regulatory_class")
            .any(|v| v.eq_ignore_ascii_case("terminator"))
    {
        Some(TerminationGeometryKind::PolymeraseTerminationRegion)
    } else if feature.kind.eq_ignore_ascii_case("polyA_signal") {
        Some(TerminationGeometryKind::ProcessingSignal)
    } else if feature.kind.eq_ignore_ascii_case("polyA_site") {
        Some(TerminationGeometryKind::AnnotatedRnaEnd)
    } else {
        None
    }
}
fn feature_kind(row: &TerminationEvidenceRow) -> &'static str {
    // Sequence-only predictions never become exact observed cleavage sites.
    match row.geometry_kind {
        TerminationGeometryKind::ProcessingSignal => "polyA_signal",
        TerminationGeometryKind::PolymeraseTerminationRegion
            if row.evidence_kind == TerminationEvidenceKind::ReferenceMatch =>
        {
            "terminator"
        }
        _ => "misc_feature",
    }
}
fn validate_catalog(catalog: &TerminationCatalog) -> Result<(), EngineError> {
    if catalog.schema != TERMINATION_CATALOG_SCHEMA
        || catalog.catalog_id.trim().is_empty()
        || catalog.version.trim().is_empty()
        || catalog.entries.is_empty()
        || catalog.entries.len() > 1000
    {
        return Err(invalid(
            "Termination catalog requires schema, identity, version and 1..1000 entries",
        ));
    }
    let mut ids = BTreeSet::new();
    let mut total_bases = 0usize;
    for entry in &catalog.entries {
        if [
            &entry.id,
            &entry.label,
            &entry.host,
            &entry.mechanism,
            &entry.source_url,
            &entry.source_version,
            &entry.recreation,
            &entry.redistribution_status,
            &entry.approval_provenance,
        ]
        .iter()
        .any(|s| s.trim().is_empty())
            || !ids.insert(&entry.id)
        {
            return Err(invalid(
                "Termination catalog entry requires unique ID and complete source/host/mechanism/redistribution provenance",
            ));
        }
        if entry.sequence.is_empty()
            || entry.sequence.len() > 100_000
            || !entry
                .sequence
                .bytes()
                .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
            || sha256_prefixed_bytes(entry.sequence.as_bytes()) != entry.sequence_sha256
            || matches!(
                entry.element_kind,
                TerminationGeometryKind::AnnotatedRnaEnd
                    | TerminationGeometryKind::CleavageWindow
                    | TerminationGeometryKind::IntrinsicTerminatorCandidate
            )
        {
            return Err(invalid(
                "Termination catalog sequence requires 1..100000 uppercase A/C/G/T bases and matching sha256: digest",
            ));
        }
        total_bases += entry.sequence.len();
        if total_bases > 1_000_000 {
            return Err(invalid(
                "Termination catalog exceeds one million sequence bases",
            ));
        }
    }
    Ok(())
}
fn finalize(mut report: TerminationScreenReport) -> Result<TerminationScreenReport, EngineError> {
    for row in &mut report.rows {
        if row.row_id.is_empty() {
            row.row_id = hash(&(&report.sequence_sha256, &report.topology, &row))?;
        }
    }
    report.rows.sort_by(|a, b| {
        (
            a.start_0based,
            a.end_0based_exclusive,
            a.reverse,
            a.evidence_kind,
            &a.element_id,
            &a.row_id,
        )
            .cmp(&(
                b.start_0based,
                b.end_0based_exclusive,
                b.reverse,
                b.evidence_kind,
                &b.element_id,
                &b.row_id,
            ))
    });
    if report.rows.len() > report.request.max_rows {
        report.rows.truncate(report.request.max_rows);
        report.truncated = true;
        report.status = "truncated".into();
        report
            .warnings
            .push("Output row cap reached; this is not a complete negative search.".into());
    }
    report.report_id = hash(&(
        &report.request,
        &report.sequence_sha256,
        &report.annotation_sha256,
        &report.catalog_sha256,
    ))?;
    report.content_sha256.clear();
    report.content_sha256 = hash(&report)?;
    Ok(report)
}

/// Retain every separately reviewed source while keeping one biological feature.
/// All parsing/serialization completes before the engine mutates any feature.
fn merge_generated_feature(
    existing: &gb_io::seq::Feature,
    incoming: gb_io::seq::Feature,
) -> Result<gb_io::seq::Feature, EngineError> {
    let mut merged = existing.clone();
    for key in [
        "termination_selected_evidence_b64",
        "termination_selected_catalog_entries_b64",
    ] {
        let mut values = BTreeMap::new();
        for feature in [existing, &incoming] {
            for text in feature.qualifier_values(key) {
                let array: Vec<serde_json::Value> = decode_provenance(text)?;
                for value in array {
                    values.insert(
                        serde_json::to_string(&value).map_err(|e| invalid(e.to_string()))?,
                        value,
                    );
                }
            }
        }
        if !values.is_empty() {
            merged.qualifiers.retain(|(k, _)| k.as_ref() != key);
            merged.qualifiers.push((
                key.into(),
                Some(encode_provenance(
                    &values.into_values().collect::<Vec<_>>(),
                )?),
            ));
        }
    }
    for (key, value) in incoming.qualifiers {
        if matches!(
            key.as_ref(),
            "label"
                | "termination_selected_evidence_b64"
                | "termination_selected_catalog_entries_b64"
        ) {
            continue;
        }
        if !merged
            .qualifiers
            .iter()
            .any(|(k, v)| *k == key && *v == value)
        {
            merged.qualifiers.push((key, value));
        }
    }
    Ok(merged)
}

impl GentleEngine {
    /// Inspect exact references, existing annotations, or explicitly selected
    /// heuristic profiles without changing sequence/project state.
    pub fn screen_transcription_termination(
        &self,
        request: &TerminationScreenRequest,
    ) -> Result<TerminationScreenReport, EngineError> {
        let dna = self
            .state
            .sequences
            .get(&request.seq_id)
            .ok_or_else(|| invalid(format!("Unknown sequence '{}'", request.seq_id)))?;
        if request.host.trim().is_empty() || request.polymerase.trim().is_empty() {
            return Err(invalid(
                "Termination screening requires explicit host and polymerase context; use unknown explicitly when genuinely unresolved",
            ));
        }
        if request.max_rows == 0
            || request.max_rows > 50_000
            || request.max_work == 0
            || request.max_work > 100_000_000
        {
            return Err(invalid(
                "Termination screen requires 1..50000 rows and 1..100000000 work units",
            ));
        }
        if request.bacterial_intrinsic
            && (!matches!(
                request.host.to_ascii_lowercase().as_str(),
                "bacteria" | "bacterial" | "e_coli" | "escherichia_coli" | "escherichia coli"
            ) || !matches!(request.polymerase.as_str(), "bacterial_rnap" | "bacterial"))
        {
            return Err(invalid(
                "Bacterial intrinsic prediction requires host=bacteria (or e_coli) and polymerase=bacterial_rnap; archaeal machinery is unsupported",
            ));
        }
        if request.mammalian_pol_ii
            && (!matches!(
                request.host.to_ascii_lowercase().as_str(),
                "mammalian" | "human" | "mouse"
            ) || !matches!(
                request.polymerase.to_ascii_lowercase().as_str(),
                "pol_ii" | "pol ii"
            ))
        {
            return Err(invalid(
                "Mammalian processing prediction requires host=mammalian (or human/mouse), polymerase=pol_ii",
            ));
        }
        let mut effective = request.clone();
        let mut strand = None;
        if let Some(id) = request.subject_feature_id {
            let feature = dna
                .features()
                .get(id)
                .ok_or_else(|| invalid("Unknown termination subject feature ID"))?;
            if !matches!(feature.kind.as_ref(), "mRNA" | "CDS" | "transcript")
                || !exact_location(&feature.location)
            {
                return Err(invalid(
                    "Termination subject must be one exact transcript/mRNA/CDS feature; imprecise or ambiguous subjects require explicit resolution",
                ));
            }
            let mut ranges = Vec::new();
            collect_location_ranges_usize(&feature.location, &mut ranges);
            let endpoint = crate::feature_location::transcript_five_prime_endpoint(feature)
                .ok_or_else(|| {
                    invalid("Termination subject has mixed/uncertain strand geometry")
                })?;
            if ranges.is_empty() || ranges.iter().any(|(s, e)| s >= e || *e > dna.len()) {
                return Err(invalid("Termination subject has out-of-sequence geometry"));
            }
            strand = Some(endpoint.reverse);
            if request.start_0based == 0 && request.end_0based_exclusive.is_none() {
                if ranges.len() != 1 {
                    return Err(invalid(
                        "Compound transcript subject requires an explicit source-DNA scan interval; mature RNA is not reconstructed",
                    ));
                }
                effective.start_0based = ranges[0].0;
                effective.end_0based_exclusive = Some(ranges[0].1);
            }
        }
        let end = effective.end_0based_exclusive.unwrap_or(dna.len());
        let start = effective.start_0based;
        if start >= end {
            return Err(invalid(
                "unsupported_span: termination scan needs ordered non-wrapping start < end",
            ));
        }
        if end > dna.len() || end - start > 2_000_000 || dna.len() > i64::MAX as usize {
            return Err(invalid(
                "Termination source-DNA interval is outside the sequence or exceeds 2 Mb",
            ));
        }
        if dna.is_circular()
            && request.end_0based_exclusive.is_none()
            && request.subject_feature_id.is_none()
        {
            return Err(invalid(
                "unsupported_span: whole-circle termination scanning needs origin-aware support; select an explicit non-wrapping interval",
            ));
        }
        effective.end_0based_exclusive = Some(end);
        let mut report = TerminationScreenReport {
            schema: TERMINATION_SCREEN_SCHEMA.into(),
            report_id: String::new(),
            content_sha256: String::new(),
            request: effective,
            sequence_sha256: sha256_prefixed_bytes(dna.forward_bytes()),
            annotation_sha256: hash(dna.features())?,
            catalog_sha256: None,
            catalog: None,
            topology: if dna.is_circular() {
                "circular"
            } else {
                "linear"
            }
            .into(),
            status: "completed".into(),
            mechanism_outcomes: vec![],
            rows: vec![],
            scanned_windows: 0,
            ambiguous_windows: 0,
            truncated: false,
            warnings: vec![],
            non_claim: TERMINATION_NON_CLAIM.into(),
        };
        if dna.is_circular() {
            report.warnings.push("Explicit non-wrapping interval on circular DNA; origin-crossing candidates were not assessed.".into());
        }
        if request.subject_feature_id.is_some() {
            report.warnings.push("Subject association is caller-selected source-DNA context, not a reconstructed mature transcript or inferred 3′ UTR.".into());
        }
        let scan = &dna.forward_bytes()[start..end];
        if request.include_annotations {
            for (id, feature) in dna.features().iter().enumerate() {
                if feature
                    .qualifier_values("gentle_generated")
                    .any(|v| v == "transcription_termination")
                {
                    continue;
                }
                let Some(geometry_kind) = annotation_geometry(feature) else {
                    continue;
                };
                let mut ranges = vec![];
                collect_location_ranges_usize(&feature.location, &mut ranges);
                if !exact_location(&feature.location) || ranges.len() != 1 {
                    report
                        .warnings
                        .push(format!("Skipped imprecise/compound source feature {id}"));
                    continue;
                }
                let (s, e) = ranges[0];
                let reverse = feature_is_reverse(feature);
                if s < start || e > end || strand.is_some_and(|v| v != reverse) {
                    continue;
                }
                if report.rows.len() > request.max_rows {
                    report.truncated = true;
                    break;
                }
                report.rows.push(TerminationEvidenceRow {
                    row_id: format!("annotation_{id}"),
                    evidence_kind: TerminationEvidenceKind::SourceAnnotation,
                    geometry_kind,
                    label: Self::feature_display_label(feature, id),
                    start_0based: s,
                    end_0based_exclusive: e,
                    reverse,
                    element_id: None,
                    source_feature_id: Some(id),
                    subject_feature_id: request.subject_feature_id,
                    sequence_5prime_to_3prime: None,
                    components: BTreeMap::new(),
                    related_intervals: vec![],
                });
            }
        }
        let mut work = 0;
        if let Some(source) = &request.catalog {
            if source.path.trim().is_empty()
                || !source.expected_sha256.starts_with("sha256:")
                || source.expected_sha256.len() != 71
            {
                return Err(invalid(
                    "Catalog path and sha256:-prefixed expected raw-byte digest are required",
                ));
            }
            match std::fs::metadata(&source.path) {
                Ok(metadata) if !metadata.is_file() => {
                    return Err(invalid("Termination catalog must be a regular file"));
                }
                _ => {}
            }
            let file = match std::fs::File::open(&source.path) {
                Ok(v) => v,
                Err(error) => {
                    report.status = "catalog_unavailable".into();
                    report.mechanism_outcomes.push(TerminationMechanismOutcome {
                        mechanism: "reference_recognition".into(),
                        status: "catalog_unavailable".into(),
                        message: format!("Cannot read declared catalog: {error}"),
                    });
                    return finalize(report);
                }
            };
            if !file
                .metadata()
                .map_err(|e| invalid(format!("Catalog metadata failure: {e}")))?
                .is_file()
            {
                return Err(invalid("Termination catalog must remain a regular file"));
            }
            let mut bytes = vec![];
            file.take((MAX_CATALOG_BYTES + 1) as u64)
                .read_to_end(&mut bytes)
                .map_err(|e| invalid(format!("Catalog read failure: {e}")))?;
            if bytes.len() > MAX_CATALOG_BYTES {
                return Err(invalid("Termination catalog exceeds 4 MiB"));
            }
            let digest = sha256_prefixed_bytes(&bytes);
            report.catalog_sha256 = Some(digest.clone());
            if digest != source.expected_sha256 {
                report.status = "freshness_refused".into();
                report.mechanism_outcomes.push(TerminationMechanismOutcome {
                    mechanism: "reference_recognition".into(),
                    status: "catalog_changed".into(),
                    message: "Catalog raw-byte digest differs from the declared expected digest"
                        .into(),
                });
                return finalize(report);
            }
            let catalog: TerminationCatalog = serde_json::from_slice(&bytes)
                .map_err(|e| invalid(format!("Malformed termination catalog: {e}")))?;
            validate_catalog(&catalog)?;
            for entry in &catalog.entries {
                if entry.host.to_ascii_lowercase() != request.host.to_ascii_lowercase()
                    && entry.host != "any"
                {
                    report.warnings.push(format!("Reference '{}' declares host '{}', requested host is '{}'; identity does not establish host compatibility",entry.id,entry.host,request.host));
                }
                let exact_scan = super::exact_reference::exact_reference_matches(
                    scan,
                    entry.sequence.as_bytes(),
                    &mut work,
                    request.max_work,
                    request.max_rows.saturating_add(1),
                    strand,
                )
                .map_err(invalid)?;
                report.scanned_windows += exact_scan.scanned_windows;
                report.ambiguous_windows += exact_scan.ambiguous_windows;
                for hit in exact_scan.hits {
                    if strand.is_some_and(|v| v != hit.reverse) {
                        continue;
                    }
                    report.rows.push(TerminationEvidenceRow {
                        row_id: format!(
                            "reference_{}_{}_{}",
                            entry.id,
                            start + hit.start,
                            hit.reverse
                        ),
                        evidence_kind: TerminationEvidenceKind::ReferenceMatch,
                        geometry_kind: entry.element_kind,
                        label: entry.label.clone(),
                        start_0based: start + hit.start,
                        end_0based_exclusive: start + hit.end,
                        reverse: hit.reverse,
                        element_id: Some(entry.id.clone()),
                        source_feature_id: None,
                        subject_feature_id: request.subject_feature_id,
                        sequence_5prime_to_3prime: Some(entry.sequence.clone()),
                        components: BTreeMap::new(),
                        related_intervals: vec![],
                    });
                    if report.rows.len() > request.max_rows {
                        report.truncated = true;
                        break;
                    }
                }
                if report.truncated {
                    break;
                }
            }
            report.catalog = Some(catalog);
            report.mechanism_outcomes.push(TerminationMechanismOutcome{mechanism:"reference_recognition".into(),status:if report.truncated {"truncated"}else{"completed"}.into(),message:"Exact sequence identity; reference provenance and host do not establish activity in this construct".into()});
        }
        if request.catalog.is_none() && !request.bacterial_intrinsic && !request.mammalian_pol_ii {
            report.status = "catalog_unavailable".into();
            report.mechanism_outcomes.push(TerminationMechanismOutcome{mechanism:"reference_recognition".into(),status:"catalog_unavailable".into(),message:"No catalog or prediction profile supplied; source annotations remain inspectable evidence".into()});
        }
        let mut missing_prediction_context = 0usize;
        if request.bacterial_intrinsic || request.mammalian_pol_ii {
            let mut prediction_request = request.clone();
            prediction_request.max_work = request.max_work.saturating_sub(work);
            let predictions =
                super::termination_predictions::predict_termination(scan, &prediction_request)
                    .map_err(invalid)?;
            missing_prediction_context = predictions.missing_context_windows;
            report.scanned_windows += predictions.scanned_windows;
            report.ambiguous_windows += predictions.ambiguous_windows;
            report.truncated |= predictions.truncated;
            report.warnings.extend(predictions.warnings);
            for mut row in predictions.rows {
                if strand.is_some_and(|v| v != row.reverse) {
                    continue;
                }
                row.start_0based += start;
                row.end_0based_exclusive += start;
                row.subject_feature_id = request.subject_feature_id;
                for related in &mut row.related_intervals {
                    related.start_0based += start;
                    related.end_0based_exclusive += start;
                }
                report.rows.push(row);
            }
        }
        for (enabled, mechanism) in [
            (request.bacterial_intrinsic, "bacterial_intrinsic"),
            (request.mammalian_pol_ii, "mammalian_pol_ii"),
        ] {
            if enabled {
                report.mechanism_outcomes.push(TerminationMechanismOutcome{mechanism:mechanism.into(),status:if report.truncated{"truncated"}else if mechanism=="mammalian_pol_ii" && missing_prediction_context>0 {"insufficient_context"}else{"completed"}.into(),message:"Authored sequence heuristic, not calibrated sensitivity/specificity or observed biological activity".into()});
            }
        }
        for mechanism in &request.requested_mechanisms {
            if report
                .mechanism_outcomes
                .iter()
                .any(|v| v.mechanism == *mechanism)
            {
                continue;
            }
            let supported = match mechanism.as_str() {
                "bacterial_intrinsic" => request.bacterial_intrinsic,
                "mammalian_pol_ii" => request.mammalian_pol_ii,
                "reference_recognition" => request.catalog.is_some(),
                _ => false,
            };
            report.mechanism_outcomes.push(TerminationMechanismOutcome {
                mechanism: mechanism.clone(),
                status: if supported {
                    "completed"
                } else if mechanism == "reference_recognition" {
                    "catalog_unavailable"
                } else {
                    "unsupported_mechanism"
                }
                .into(),
                message: if supported {
                    "Requested analysis mode enabled"
                } else if mechanism == "reference_recognition" {
                    "No catalog supplied for requested reference recognition"
                } else {
                    "No de novo detector is enabled for this mechanism; no absence claim follows"
                }
                .into(),
            });
        }
        if report.truncated {
            report.status = "truncated".into();
        }
        finalize(report)
    }

    /// Recompute and verify the reviewed report before one atomic feature append.
    pub fn materialize_transcription_termination(
        &mut self,
        request: &TerminationMaterializeRequest,
    ) -> Result<TerminationScreenReport, EngineError> {
        let report = self.screen_transcription_termination(&request.screen)?;
        if report.status != "completed" || report.truncated {
            return Err(invalid(
                "Cannot annotate an unavailable, refused or incomplete termination screen",
            ));
        }
        if report.content_sha256 != request.expected_report_sha256 {
            return Err(invalid(
                "Termination evidence changed since review; inspect again before adding features",
            ));
        }
        let selected = request.row_ids.iter().collect::<BTreeSet<_>>();
        if selected.is_empty()
            || selected.len() != request.row_ids.len()
            || selected
                .iter()
                .any(|id| !report.rows.iter().any(|r| &r.row_id == *id))
        {
            return Err(invalid(
                "Select unique existing termination evidence row IDs",
            ));
        }
        let dna = self
            .state
            .sequences
            .get(&request.screen.seq_id)
            .ok_or_else(|| invalid("Termination source missing"))?;
        let mut additions = Vec::new();
        let mut updates = Vec::new();
        let mut groups: BTreeMap<(&str, usize, usize, bool), Vec<&TerminationEvidenceRow>> =
            BTreeMap::new();
        for row in report.rows.iter().filter(|r| {
            selected.contains(&r.row_id)
                && r.evidence_kind != TerminationEvidenceKind::SourceAnnotation
        }) {
            groups
                .entry((
                    feature_kind(row),
                    row.start_0based,
                    row.end_0based_exclusive,
                    row.reverse,
                ))
                .or_default()
                .push(row);
        }
        for ((kind, _, _, _), rows) in groups {
            let row = rows[0];
            let keys = rows
                .iter()
                .map(|row| {
                    hash(&(
                        &report.sequence_sha256,
                        &report.topology,
                        row,
                        &report.catalog_sha256,
                    ))
                })
                .collect::<Result<Vec<_>, _>>()?;
            let key = keys[0].clone();
            if dna.features().iter().any(|feature| {
                keys.iter().all(|key| {
                    feature
                        .qualifier_values("termination_evidence_id")
                        .any(|v| compact_digest(v) == *key)
                })
            }) {
                continue;
            }
            let same_geometry = |feature: &gb_io::seq::Feature| {
                let equivalent = feature.kind.eq_ignore_ascii_case(kind)
                    || (kind == "terminator"
                        && annotation_geometry(feature)
                            == Some(TerminationGeometryKind::PolymeraseTerminationRegion));
                let mut ranges = vec![];
                collect_location_ranges_usize(&feature.location, &mut ranges);
                equivalent
                    && exact_location(&feature.location)
                    && ranges == vec![(row.start_0based, row.end_0based_exclusive)]
                    && feature_is_reverse(feature) == row.reverse
            };
            let is_generated = |feature: &gb_io::seq::Feature| {
                feature
                    .qualifier_values("gentle_generated")
                    .any(|v| v == "transcription_termination")
            };
            if kind != "misc_feature"
                && dna
                    .features()
                    .iter()
                    .any(|feature| same_geometry(feature) && !is_generated(feature))
            {
                continue;
            }
            let existing_generated = dna
                .features()
                .iter()
                .position(|feature| same_geometry(feature) && is_generated(feature));
            let mut location = gb_io::seq::Location::simple_range(
                row.start_0based as i64,
                row.end_0based_exclusive as i64,
            );
            if row.reverse {
                location = gb_io::seq::Location::Complement(Box::new(location));
            }
            let mut qualifiers = vec![
                ("label".into(), Some(format!("Candidate {}", row.label))),
                (
                    "gentle_generated".into(),
                    Some("transcription_termination".into()),
                ),
                ("termination_evidence_id".into(), Some(key)),
                (
                    "termination_selected_evidence_b64".into(),
                    Some(encode_provenance(&rows)?),
                ),
                (
                    "termination_report_sha256".into(),
                    Some(report.content_sha256.clone()),
                ),
                (
                    "termination_evidence_b64".into(),
                    Some(encode_provenance(row)?),
                ),
                (
                    "termination_request_b64".into(),
                    Some(encode_provenance(&report.request)?),
                ),
                ("note".into(), Some(TERMINATION_NON_CLAIM.into())),
            ];
            for key in keys.iter().skip(1) {
                qualifiers.push(("termination_evidence_id".into(), Some(key.clone())));
            }
            let selected_entries =
                rows.iter()
                    .filter_map(|row| row.element_id.as_ref())
                    .filter_map(|id| {
                        report.catalog.as_ref().and_then(|catalog| {
                            catalog.entries.iter().find(|entry| &entry.id == id)
                        })
                    })
                    .collect::<Vec<_>>();
            if !selected_entries.is_empty() {
                qualifiers.push((
                    "termination_selected_catalog_entries_b64".into(),
                    Some(encode_provenance(&selected_entries)?),
                ));
                for entry in selected_entries {
                    qualifiers.push((
                        "termination_catalog_entry_b64".into(),
                        Some(encode_provenance(entry)?),
                    ));
                }
                qualifiers.push((
                    "termination_catalog_sha256".into(),
                    report.catalog_sha256.clone(),
                ));
            }
            if kind == "terminator" {
                qualifiers.push(("regulatory_class".into(), Some("terminator".into())));
            }
            let incoming = gb_io::seq::Feature {
                kind: kind.into(),
                location,
                qualifiers,
            };
            if let Some(index) = existing_generated {
                updates.push((
                    index,
                    merge_generated_feature(&dna.features()[index], incoming)?,
                ));
            } else {
                additions.push(incoming);
            }
        }
        let dna = self
            .state
            .sequences
            .get_mut(&request.screen.seq_id)
            .ok_or_else(|| invalid("Termination source missing"))?;
        if !additions.is_empty() || !updates.is_empty() {
            for (index, feature) in updates {
                dna.features_mut()[index] = feature;
            }
            dna.features_mut().extend(additions);
            Self::prepare_sequence(dna);
        }
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gb_io::seq::{Feature, Location};
    // Synthetic exact element, source annotations and DNA authored here. They
    // are recreated by these tests, used only for software/provenance boundaries,
    // and do not represent a validated biological terminator or licensed asset.
    fn setup() -> (GentleEngine, tempfile::TempDir, TerminationScreenRequest) {
        let mut engine = GentleEngine::new();
        engine.state_mut().sequences.insert(
            "synthetic".into(),
            DNAsequence::from_sequence("CCAACGTTGGCCAACGTT").unwrap(),
        );
        let dir = tempfile::tempdir().unwrap();
        let entry = TerminationCatalogEntry {
            id: "toy".into(),
            label: "synthetic element".into(),
            sequence: "AACG".into(),
            sequence_sha256: sha256_prefixed_bytes(b"AACG"),
            host: "bacteria".into(),
            mechanism: "bacterial_intrinsic".into(),
            element_kind: TerminationGeometryKind::PolymeraseTerminationRegion,
            source_url: "synthetic:src/engine/analysis/transcription_termination.rs".into(),
            source_version: "1".into(),
            recreation: "Run inline setup; authored canonical sequence".into(),
            redistribution_status: "synthetic".into(),
            approval_provenance: "Synthetic test only; no biological suitability assertion".into(),
        };
        let bytes = serde_json::to_vec(&TerminationCatalog {
            schema: TERMINATION_CATALOG_SCHEMA.into(),
            catalog_id: "toy".into(),
            version: "1".into(),
            entries: vec![entry],
        })
        .unwrap();
        let path = dir.path().join("catalog with spaces.json");
        std::fs::write(&path, &bytes).unwrap();
        let request = TerminationScreenRequest {
            seq_id: "synthetic".into(),
            host: "bacteria".into(),
            polymerase: "bacterial".into(),
            catalog: Some(TerminationCatalogSource {
                path: path.to_string_lossy().into(),
                expected_sha256: sha256_prefixed_bytes(&bytes),
            }),
            ..Default::default()
        };
        (engine, dir, request)
    }
    fn materialize(report: &TerminationScreenReport) -> TerminationMaterializeRequest {
        TerminationMaterializeRequest {
            screen: report.request.clone(),
            expected_report_sha256: report.content_sha256.clone(),
            row_ids: report
                .rows
                .iter()
                .filter(|r| r.evidence_kind != TerminationEvidenceKind::SourceAnnotation)
                .map(|r| r.row_id.clone())
                .collect(),
        }
    }
    #[test]
    fn termination_exact_read_only_hashes_strands_and_catalog_provenance() {
        let (engine, _dir, request) = setup();
        let before = serde_json::to_value(engine.state()).unwrap();
        let report = engine.screen_transcription_termination(&request).unwrap();
        assert_eq!(
            report,
            engine.screen_transcription_termination(&request).unwrap()
        );
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
        assert!(report.rows.iter().any(|r| r.reverse));
        assert!(report.rows.iter().any(|r| !r.reverse));
        assert!(report.catalog.is_some());
        let mut empty = report.clone();
        empty.content_sha256.clear();
        assert_eq!(report.content_sha256, hash(&empty).unwrap());
        assert_eq!(
            serde_json::from_str::<TerminationScreenReport>(
                &serde_json::to_string(&report).unwrap()
            )
            .unwrap(),
            report
        );
    }
    #[test]
    fn termination_stale_catalog_absent_malformed_and_zero_hits_are_distinct() {
        let (engine, _dir, mut request) = setup();
        let report = engine.screen_transcription_termination(&request).unwrap();
        let source = request.catalog.as_mut().unwrap();
        std::fs::write(&source.path, b"{}").unwrap();
        assert_eq!(
            engine
                .screen_transcription_termination(&request)
                .unwrap()
                .status,
            "freshness_refused"
        );
        let mut zero_request = report.request.clone();
        zero_request.start_0based = 0;
        zero_request.end_0based_exclusive = Some(2);
        // Restore the bound valid catalog before testing a genuine completed zero-hit interval.
        let valid = serde_json::to_vec(report.catalog.as_ref().unwrap()).unwrap();
        std::fs::write(&zero_request.catalog.as_ref().unwrap().path, &valid).unwrap();
        zero_request.catalog.as_mut().unwrap().expected_sha256 = sha256_prefixed_bytes(&valid);
        let zero = engine
            .screen_transcription_termination(&zero_request)
            .unwrap();
        assert_eq!(zero.status, "completed");
        assert!(zero.rows.is_empty());
        std::fs::write(&zero_request.catalog.as_ref().unwrap().path, b"{}").unwrap();
        request.catalog.as_mut().unwrap().expected_sha256 = sha256_prefixed_bytes(b"{}");
        assert!(engine.screen_transcription_termination(&request).is_err());
        std::fs::remove_file(&request.catalog.as_ref().unwrap().path).unwrap();
        assert_eq!(
            engine
                .screen_transcription_termination(&request)
                .unwrap()
                .status,
            "catalog_unavailable"
        );
        let mut source_only = report.request;
        source_only.catalog = None;
        source_only.requested_mechanisms = vec!["rho_dependent".into()];
        let r = engine
            .screen_transcription_termination(&source_only)
            .unwrap();
        assert_eq!(r.status, "catalog_unavailable");
        assert!(r.rows.is_empty());
        assert!(
            r.mechanism_outcomes
                .iter()
                .any(|v| v.status == "unsupported_mechanism")
        );
    }
    #[test]
    fn termination_apply_is_review_bound_idempotent_and_undoable() {
        let (mut engine, _dir, request) = setup();
        let report = engine.screen_transcription_termination(&request).unwrap();
        let before = engine.state().sequences["synthetic"].features().len();
        let mut bad = materialize(&report);
        bad.row_ids.push("missing".into());
        assert!(engine.materialize_transcription_termination(&bad).is_err());
        assert_eq!(
            before,
            engine.state().sequences["synthetic"].features().len()
        );
        let result = engine
            .apply(Operation::MaterializeTranscriptionTerminationFeatures {
                request: materialize(&report),
            })
            .unwrap();
        assert!(result.changed_seq_ids.contains(&"synthetic".into()));
        let after = engine.state().sequences["synthetic"].features().len();
        assert!(after > before);
        assert!(
            engine
                .materialize_transcription_termination(&materialize(&report))
                .is_err()
        ); // annotation identity changed
        let fresh = engine.screen_transcription_termination(&request).unwrap();
        engine
            .materialize_transcription_termination(&materialize(&fresh))
            .unwrap();
        assert_eq!(
            after,
            engine.state().sequences["synthetic"].features().len()
        );
        engine.undo_last_operation().unwrap();
        assert_eq!(
            before,
            engine.state().sequences["synthetic"].features().len()
        );
        engine.redo_last_operation().unwrap();
        assert_eq!(
            after,
            engine.state().sequences["synthetic"].features().len()
        );
    }
    #[test]
    fn termination_existing_regulatory_class_is_preserved_without_duplicate() {
        let (mut engine, _dir, request) = setup();
        engine
            .state_mut()
            .sequences
            .get_mut("synthetic")
            .unwrap()
            .features_mut()
            .push(Feature {
                kind: "terminator".into(),
                location: Location::simple_range(2, 6),
                qualifiers: vec![("regulatory_class".into(), Some("terminator".into()))],
            });
        let report = engine.screen_transcription_termination(&request).unwrap();
        assert!(
            report
                .rows
                .iter()
                .any(|r| r.evidence_kind == TerminationEvidenceKind::SourceAnnotation)
        );
        engine
            .materialize_transcription_termination(&materialize(&report))
            .unwrap();
        assert_eq!(
            engine.state().sequences["synthetic"]
                .features()
                .iter()
                .filter(|f| {
                    let mut ranges = vec![];
                    collect_location_ranges_usize(&f.location, &mut ranges);
                    ranges == vec![(2, 6)] && !feature_is_reverse(f)
                })
                .count(),
            1
        );
        let exported = engine.state().sequences["synthetic"]
            .to_genbank_string()
            .unwrap();
        assert!(exported.contains("termination_catalog_entry_b64"));
        let seq = gb_io::reader::SeqReader::new(exported.as_bytes())
            .next()
            .unwrap()
            .unwrap();
        let imported = DNAsequence::from_genbank_seq(seq);
        let before = imported.features().len();
        engine
            .state_mut()
            .sequences
            .insert("synthetic".into(), imported);
        let fresh = engine.screen_transcription_termination(&request).unwrap();
        engine
            .materialize_transcription_termination(&materialize(&fresh))
            .unwrap();
        assert_eq!(
            engine.state().sequences["synthetic"].features().len(),
            before
        );
    }
    #[test]
    fn termination_colocated_reference_provenance_is_aggregated() {
        let (mut engine, _dir, mut request) = setup();
        let source = request.catalog.as_mut().unwrap();
        let mut catalog: TerminationCatalog =
            serde_json::from_slice(&std::fs::read(&source.path).unwrap()).unwrap();
        let mut second = catalog.entries[0].clone();
        second.id = "second_source".into();
        second.source_url = "synthetic:second independently declared source".into();
        second.recreation =
            "Preserve  multiple spaces, \"quotes\", β and\nnewlines in source provenance".into();
        catalog.entries.push(second);
        let bytes = serde_json::to_vec(&catalog).unwrap();
        std::fs::write(&source.path, &bytes).unwrap();
        source.expected_sha256 = sha256_prefixed_bytes(&bytes);
        let report = engine.screen_transcription_termination(&request).unwrap();
        engine
            .materialize_transcription_termination(&materialize(&report))
            .unwrap();
        for feature in engine.state().sequences["synthetic"].features() {
            let entries = feature
                .qualifier_values("termination_selected_catalog_entries_b64")
                .next()
                .unwrap();
            let entries: Vec<TerminationCatalogEntry> = decode_provenance(entries).unwrap();
            assert_eq!(entries.len(), 2);
            assert!(entries.iter().any(|entry| entry.id == "second_source"));
            let rows = feature
                .qualifier_values("termination_selected_evidence_b64")
                .next()
                .unwrap();
            assert_eq!(
                decode_provenance::<Vec<TerminationEvidenceRow>>(rows)
                    .unwrap()
                    .len(),
                2
            );
        }
        let exported = engine.state().sequences["synthetic"]
            .to_genbank_string()
            .unwrap();
        let parsed = gb_io::reader::SeqReader::new(exported.as_bytes())
            .next()
            .unwrap()
            .unwrap();
        for feature in &parsed.features {
            let text = feature
                .qualifier_values("termination_selected_catalog_entries_b64")
                .next()
                .unwrap();
            let entries: Vec<TerminationCatalogEntry> = decode_provenance(text).unwrap();
            assert!(entries.iter().any(|entry| entry.id == "second_source"));
            for entry in &entries {
                assert_eq!(
                    entry,
                    catalog
                        .entries
                        .iter()
                        .find(|expected| expected.id == entry.id)
                        .unwrap()
                );
            }
            let row_text = feature
                .qualifier_values("termination_selected_evidence_b64")
                .next()
                .unwrap();
            let decoded_rows: Vec<TerminationEvidenceRow> = decode_provenance(row_text).unwrap();
            for row in &decoded_rows {
                assert_eq!(
                    row,
                    report
                        .rows
                        .iter()
                        .find(|expected| expected.row_id == row.row_id)
                        .unwrap()
                );
            }
        }
    }
    #[test]
    fn termination_sequential_colocated_reviews_merge_and_undo_provenance() {
        let (mut engine, _dir, mut request) = setup();
        let source = request.catalog.as_mut().unwrap();
        let mut catalog: TerminationCatalog =
            serde_json::from_slice(&std::fs::read(&source.path).unwrap()).unwrap();
        let mut second = catalog.entries[0].clone();
        second.id = "sequential_second".into();
        second.source_url = "synthetic:second reviewed source".into();
        catalog.entries.push(second);
        let bytes = serde_json::to_vec(&catalog).unwrap();
        std::fs::write(&source.path, &bytes).unwrap();
        source.expected_sha256 = sha256_prefixed_bytes(&bytes);
        let first = engine.screen_transcription_termination(&request).unwrap();
        let mut first_apply = materialize(&first);
        first_apply.row_ids.retain(|id| {
            first
                .rows
                .iter()
                .any(|row| &row.row_id == id && row.element_id.as_deref() == Some("toy"))
        });
        engine
            .apply(Operation::MaterializeTranscriptionTerminationFeatures {
                request: first_apply,
            })
            .unwrap();
        let count = engine.state().sequences["synthetic"].features().len();
        let before = engine.state().sequences["synthetic"].features().to_vec();
        let second = engine.screen_transcription_termination(&request).unwrap();
        let mut second_apply = materialize(&second);
        second_apply.row_ids.retain(|id| {
            second.rows.iter().any(|row| {
                &row.row_id == id && row.element_id.as_deref() == Some("sequential_second")
            })
        });
        let applied = engine
            .apply(Operation::MaterializeTranscriptionTerminationFeatures {
                request: second_apply,
            })
            .unwrap();
        assert!(applied.changed_seq_ids.contains(&"synthetic".into()));
        assert_eq!(
            engine.state().sequences["synthetic"].features().len(),
            count
        );
        for feature in engine.state().sequences["synthetic"].features() {
            assert_eq!(
                feature.qualifier_values("termination_evidence_id").count(),
                2
            );
            let entries = feature
                .qualifier_values("termination_selected_catalog_entries_b64")
                .next()
                .unwrap();
            assert_eq!(
                decode_provenance::<Vec<TerminationCatalogEntry>>(entries)
                    .unwrap()
                    .len(),
                2
            );
        }
        engine.undo_last_operation().unwrap();
        assert_eq!(
            serde_json::to_value(engine.state().sequences["synthetic"].features()).unwrap(),
            serde_json::to_value(before).unwrap()
        );
    }
    #[test]
    fn termination_catalog_rejects_nonregular_directory_without_reading() {
        let (engine, dir, mut request) = setup();
        request.catalog.as_mut().unwrap().path = dir.path().to_string_lossy().into();
        let error = engine
            .screen_transcription_termination(&request)
            .unwrap_err();
        assert!(error.message.contains("regular file"));
    }
    #[test]
    fn termination_subject_binding_source_intervals_and_origin_refusal() {
        let (mut engine, _dir, mut request) = setup();
        engine
            .state_mut()
            .sequences
            .get_mut("synthetic")
            .unwrap()
            .features_mut()
            .push(Feature {
                kind: "CDS".into(),
                location: Location::Complement(Box::new(Location::simple_range(10, 18))),
                qualifiers: vec![],
            });
        request.subject_feature_id = Some(0);
        request.end_0based_exclusive = Some(10);
        let report = engine.screen_transcription_termination(&request).unwrap();
        assert!(
            report
                .rows
                .iter()
                .all(|r| r.reverse && r.subject_feature_id == Some(0))
        );
        // Explicit scan intervals do not legitimize ambiguous or out-of-bounds
        // subject annotations. These are synthetic feature geometries only.
        engine
            .state_mut()
            .sequences
            .get_mut("synthetic")
            .unwrap()
            .features_mut()[0]
            .location = Location::Join(vec![
            Location::simple_range(2, 6),
            Location::Complement(Box::new(Location::simple_range(8, 12))),
        ]);
        assert!(engine.screen_transcription_termination(&request).is_err());
        engine
            .state_mut()
            .sequences
            .get_mut("synthetic")
            .unwrap()
            .features_mut()[0]
            .location = Location::simple_range(2, 999);
        assert!(engine.screen_transcription_termination(&request).is_err());
        request.subject_feature_id = None;
        request.start_0based = 12;
        request.end_0based_exclusive = Some(4);
        assert!(engine.screen_transcription_termination(&request).is_err());
        request.start_0based = 0;
        request.end_0based_exclusive = None;
        engine
            .state_mut()
            .sequences
            .get_mut("synthetic")
            .unwrap()
            .set_circular(true);
        assert!(engine.screen_transcription_termination(&request).is_err());
        request.end_0based_exclusive = Some(10);
        assert!(
            engine
                .screen_transcription_termination(&request)
                .unwrap()
                .warnings
                .iter()
                .any(|w| w.contains("origin-crossing"))
        );
    }
    #[test]
    fn termination_catalog_raw_lf_crlf_bytes_are_not_normalized() {
        let (engine, _dir, mut request) = setup();
        let source = request.catalog.as_mut().unwrap();
        let catalog: TerminationCatalog =
            serde_json::from_slice(&std::fs::read(&source.path).unwrap()).unwrap();
        let lf = serde_json::to_string_pretty(&catalog).unwrap() + "\n";
        std::fs::write(&source.path, lf.as_bytes()).unwrap();
        source.expected_sha256 = sha256_prefixed_bytes(lf.as_bytes());
        let lf_report = engine.screen_transcription_termination(&request).unwrap();
        let crlf = lf.replace('\n', "\r\n");
        std::fs::write(&request.catalog.as_ref().unwrap().path, crlf.as_bytes()).unwrap();
        assert_eq!(
            engine
                .screen_transcription_termination(&request)
                .unwrap()
                .status,
            "freshness_refused"
        );
        request.catalog.as_mut().unwrap().expected_sha256 = sha256_prefixed_bytes(crlf.as_bytes());
        let crlf_report = engine.screen_transcription_termination(&request).unwrap();
        assert_eq!(crlf_report.status, "completed");
        assert_ne!(lf_report.catalog_sha256, crlf_report.catalog_sha256);
        assert_eq!(lf_report.rows, crlf_report.rows);
    }

    #[test]
    fn termination_digest_does_not_depend_on_rnafold_resolution() {
        let (engine, dir, mut request) = setup();
        request.bacterial_intrinsic = true;
        let missing = dir.path().join("missing-rnafold");
        let present = dir.path().join("declared-rnafold");
        // Synthetic inert local file, not an installed folding program. The
        // detector must not inspect or execute either configured tool path.
        std::fs::write(&present, b"not a folding executable").unwrap();
        let missing_digest = {
            let _tool = crate::tool_overrides::ScopedToolOverrideGuard::set(
                crate::rna_structure::RNAFOLD_ENV_BIN,
                &missing.to_string_lossy(),
            );
            engine
                .screen_transcription_termination(&request)
                .unwrap()
                .content_sha256
        };
        let present_digest = {
            let _tool = crate::tool_overrides::ScopedToolOverrideGuard::set(
                crate::rna_structure::RNAFOLD_ENV_BIN,
                &present.to_string_lossy(),
            );
            engine
                .screen_transcription_termination(&request)
                .unwrap()
                .content_sha256
        };
        assert_eq!(missing_digest, present_digest);
        assert_eq!(std::fs::read(present).unwrap(), b"not a folding executable");
    }

    #[test]
    fn termination_work_and_row_caps_do_not_claim_absence() {
        let (engine, _dir, mut request) = setup();
        request.max_work = 1;
        assert!(engine.screen_transcription_termination(&request).is_err());
        request.max_work = 10000;
        request.max_rows = 1;
        let report = engine.screen_transcription_termination(&request).unwrap();
        assert!(report.truncated);
        assert_eq!(report.status, "truncated");
        let mut engine = engine;
        assert!(
            engine
                .materialize_transcription_termination(&materialize(&report))
                .is_err()
        );
    }
}
