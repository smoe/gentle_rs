//! Reporter catalog validation, recommendation, and corpus export operations.

use super::*;
use sha1::{Digest, Sha1};

const DEFAULT_REPORTER_CATALOG_PATH: &str = "assets/reporter_catalog.json";
const DEFAULT_REPORTER_RECOMMENDATION_LIMIT: usize = 5;

impl GentleEngine {
    pub(crate) fn list_reporter_catalog(
        &self,
        catalog_path: Option<&str>,
        filter: Option<&str>,
        limit: Option<usize>,
    ) -> Result<ReporterCatalogReport, EngineError> {
        let mut report = self.annotated_reporter_catalog(catalog_path, &[])?;
        report.records = Self::filter_reporter_records(report.records, filter, limit);
        report.active_record_count = report.records.len();
        report.record_count = report.active_record_count + report.quarantined_record_count;
        Ok(report)
    }

    pub(crate) fn recommend_reporters(
        &self,
        catalog_path: Option<&str>,
        constraints: ReporterConstraints,
        limit: Option<usize>,
    ) -> Result<ReporterRecommendationResult, EngineError> {
        let normalized_constraints = Self::normalize_reporter_constraints(constraints)?;
        let catalog_report = self
            .annotated_reporter_catalog(catalog_path, &normalized_constraints.forbidden_motifs)?;
        let mut accepted = Vec::new();
        let mut rejected = Vec::new();

        for entry in catalog_report.records {
            let reasons = Self::reporter_hard_rejection_reasons(&entry, &normalized_constraints)?;
            if reasons.is_empty() {
                accepted.push(Self::score_reporter_candidate(
                    &entry,
                    &normalized_constraints,
                ));
            } else {
                rejected.push(ReporterRejectedCandidate {
                    reporter_id: entry.record.id,
                    name: entry.record.name,
                    reasons,
                });
            }
        }

        accepted.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(Ordering::Equal)
                .then_with(|| a.name.cmp(&b.name))
                .then_with(|| a.reporter_id.cmp(&b.reporter_id))
        });
        let accepted_candidate_count = accepted.len();
        let requested_limit = limit
            .unwrap_or(DEFAULT_REPORTER_RECOMMENDATION_LIMIT)
            .max(1);
        accepted.truncate(requested_limit);
        for (idx, row) in accepted.iter_mut().enumerate() {
            row.rank = idx + 1;
        }

        let mut warnings = catalog_report.warnings;
        if !catalog_report.quarantined_records.is_empty() {
            warnings.push(format!(
                "{} reporter catalog row(s) were quarantined before ranking",
                catalog_report.quarantined_records.len()
            ));
        }
        if accepted.is_empty() {
            warnings.push(
                "No reporter candidate satisfied all hard constraints; inspect rejected_candidates"
                    .to_string(),
            );
        }

        Ok(ReporterRecommendationResult {
            schema: REPORTER_RECOMMENDATION_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            catalog_path: Self::effective_reporter_catalog_label(catalog_path),
            constraints: normalized_constraints,
            considered_candidate_count: accepted_candidate_count + rejected.len(),
            recommended_candidate_count: accepted.len(),
            rejected_candidate_count: rejected.len(),
            recommendations: accepted,
            rejected_candidates: rejected,
            warnings,
        })
    }

    pub(crate) fn export_reporter_corpus(
        &self,
        catalog_path: Option<&str>,
        path: &str,
        format: ReporterCorpusExportFormat,
    ) -> Result<ReporterCorpusExport, EngineError> {
        let catalog = self.annotated_reporter_catalog(catalog_path, &[])?;
        let export = ReporterCorpusExport {
            schema: REPORTER_CORPUS_EXPORT_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            catalog_path: catalog.catalog_path,
            format,
            record_count: catalog.records.len(),
            records: catalog.records,
            warnings: catalog.warnings,
        };
        match format {
            ReporterCorpusExportFormat::Json => {
                self.write_pretty_json_file(&export, path, "reporter corpus export")?;
            }
            ReporterCorpusExportFormat::Jsonl => {
                let mut lines = String::new();
                for record in &export.records {
                    let line = serde_json::to_string(record).map_err(|e| EngineError {
                        code: ErrorCode::Internal,
                        message: format!("Could not serialize reporter corpus JSONL row: {e}"),
                        cause_chain: vec![],
                    })?;
                    lines.push_str(&line);
                    lines.push('\n');
                }
                self.write_text_file(path, &lines, "reporter corpus JSONL export")?;
            }
        }
        Ok(export)
    }

    fn annotated_reporter_catalog(
        &self,
        catalog_path: Option<&str>,
        forbidden_motifs: &[String],
    ) -> Result<ReporterCatalogReport, EngineError> {
        let catalog_label = Self::effective_reporter_catalog_label(catalog_path);
        let catalog = Self::load_reporter_catalog(catalog_path)?;
        if catalog.schema != REPORTER_CATALOG_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Reporter catalog '{}' has schema '{}', expected '{}'",
                    catalog_label, catalog.schema, REPORTER_CATALOG_SCHEMA
                ),
                cause_chain: vec![],
            });
        }

        let mut seen_ids = BTreeSet::new();
        let mut records = Vec::new();
        let mut quarantined_records = Vec::new();
        let mut warnings = Vec::new();
        for mut record in catalog.records {
            record.sequence = Self::normalize_reporter_sequence(&record.sequence);
            let mut reasons = Self::reporter_catalog_quarantine_reasons(&record);
            if !record.id.trim().is_empty() && !seen_ids.insert(record.id.clone()) {
                reasons.push("duplicate_reporter_id".to_string());
            }
            let annotation = Self::annotate_reporter_record(&record, forbidden_motifs)?;
            if !annotation.checksum_ok {
                reasons.push("sequence_sha1_mismatch".to_string());
            }
            if reasons.is_empty() {
                let mut row_warnings = Vec::new();
                if !annotation.likely_complete_cds {
                    row_warnings.push(
                        "sequence is not a complete CDS by simple ATG/frame/stop checks"
                            .to_string(),
                    );
                }
                records.push(ReporterAnnotatedRecord {
                    record,
                    annotation,
                    warnings: row_warnings,
                });
            } else {
                quarantined_records.push(ReporterQuarantinedRecord {
                    id: record.id,
                    name: record.name,
                    reasons,
                });
            }
        }
        records.sort_by(|a, b| a.record.id.cmp(&b.record.id));
        quarantined_records.sort_by(|a, b| a.id.cmp(&b.id));
        if !quarantined_records.is_empty() {
            warnings.push(format!(
                "{} reporter catalog row(s) quarantined by provenance/safety gates",
                quarantined_records.len()
            ));
        }
        let active_record_count = records.len();
        Ok(ReporterCatalogReport {
            schema: REPORTER_CATALOG_REPORT_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            catalog_path: catalog_label,
            record_count: active_record_count + quarantined_records.len(),
            active_record_count,
            quarantined_record_count: quarantined_records.len(),
            records,
            quarantined_records,
            warnings,
        })
    }

    fn load_reporter_catalog(catalog_path: Option<&str>) -> Result<ReporterCatalog, EngineError> {
        let label = Self::effective_reporter_catalog_label(catalog_path);
        let text = match catalog_path.map(str::trim).filter(|path| !path.is_empty()) {
            None => include_str!("../../../assets/reporter_catalog.json").to_string(),
            Some(path) if path == DEFAULT_REPORTER_CATALOG_PATH => {
                include_str!("../../../assets/reporter_catalog.json").to_string()
            }
            Some(path) => std::fs::read_to_string(path).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not read reporter catalog '{path}': {e}"),
                cause_chain: vec![],
            })?,
        };
        serde_json::from_str(&text).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not parse reporter catalog '{}': {e}", label),
            cause_chain: vec![],
        })
    }

    fn effective_reporter_catalog_label(catalog_path: Option<&str>) -> String {
        catalog_path
            .map(str::trim)
            .filter(|path| !path.is_empty())
            .unwrap_or(DEFAULT_REPORTER_CATALOG_PATH)
            .to_string()
    }

    fn reporter_catalog_quarantine_reasons(record: &ReporterRecord) -> Vec<String> {
        let mut reasons = Vec::new();
        if record.id.trim().is_empty() {
            reasons.push("missing_id".to_string());
        }
        if record.name.trim().is_empty() {
            reasons.push("missing_name".to_string());
        }
        if record.sequence.trim().is_empty() {
            reasons.push("missing_sequence".to_string());
        } else if !Self::is_valid_reporter_sequence(&record.sequence) {
            reasons.push("invalid_sequence_characters".to_string());
        }
        if record.sequence_sha1.trim().is_empty() {
            reasons.push("missing_sequence_sha1".to_string());
        }
        if !record.license_status.starts_with("open") {
            reasons.push("unclear_or_unaccepted_license_status".to_string());
        }
        if record.source_refs.is_empty() {
            reasons.push("missing_source_reference".to_string());
        }
        if record.provenance_note.trim().is_empty() {
            reasons.push("missing_provenance_note".to_string());
        }
        if record.safety_scope != "benign_reporter_only" {
            reasons.push("outside_benign_reporter_scope".to_string());
        }
        if !matches!(
            record.reporter_class.as_str(),
            "fluorescent_protein" | "chromoprotein" | "luciferase" | "colorimetric_enzyme"
        ) {
            reasons.push("unsupported_reporter_class".to_string());
        }
        if record
            .source_refs
            .iter()
            .any(|source| !source.license_status.starts_with("open"))
        {
            reasons.push("source_reference_license_not_open".to_string());
        }
        reasons
    }

    fn annotate_reporter_record(
        record: &ReporterRecord,
        forbidden_motifs: &[String],
    ) -> Result<ReporterComputedAnnotation, EngineError> {
        let sequence = record.sequence.as_bytes();
        let gc_fraction = Self::sequence_gc_fraction(sequence);
        let starts_with_atg = sequence
            .get(0..3)
            .map(|start| start.eq_ignore_ascii_case(b"ATG"))
            .unwrap_or(false);
        let ends_with_stop = sequence
            .get(sequence.len().saturating_sub(3)..sequence.len())
            .map(|stop| matches!(stop, b"TAA" | b"TAG" | b"TGA"))
            .unwrap_or(false);
        let multiple_of_three = sequence.len() % 3 == 0;
        let checksum = format!("{:x}", Sha1::digest(record.sequence.as_bytes()));
        let mut forbidden_motif_hits = Vec::new();
        for motif in forbidden_motifs {
            if Self::contains_motif_any_strand(sequence, motif)? {
                forbidden_motif_hits.push(motif.clone());
            }
        }
        Ok(ReporterComputedAnnotation {
            length_bp: sequence.len(),
            gc_fraction,
            starts_with_atg,
            ends_with_stop,
            multiple_of_three,
            likely_complete_cds: starts_with_atg && multiple_of_three && ends_with_stop,
            checksum_ok: checksum == record.sequence_sha1.to_ascii_lowercase(),
            forbidden_motif_hits,
        })
    }

    fn normalize_reporter_sequence(sequence: &str) -> String {
        sequence
            .as_bytes()
            .iter()
            .filter(|b| !b.is_ascii_whitespace())
            .map(|b| match b.to_ascii_uppercase() {
                b'U' => 'T',
                other => other as char,
            })
            .collect()
    }

    fn is_valid_reporter_sequence(sequence: &str) -> bool {
        sequence
            .as_bytes()
            .iter()
            .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
    }

    fn filter_reporter_records(
        records: Vec<ReporterAnnotatedRecord>,
        filter: Option<&str>,
        limit: Option<usize>,
    ) -> Vec<ReporterAnnotatedRecord> {
        let filter = filter
            .map(Self::normalize_reporter_token)
            .filter(|value| !value.is_empty());
        let mut out = records
            .into_iter()
            .filter(|row| match filter.as_deref() {
                None => true,
                Some(token) => Self::reporter_search_text(row).contains(token),
            })
            .collect::<Vec<_>>();
        if let Some(limit) = limit {
            out.truncate(limit);
        }
        out
    }

    fn reporter_search_text(row: &ReporterAnnotatedRecord) -> String {
        let mut parts = vec![
            row.record.id.clone(),
            row.record.name.clone(),
            row.record.reporter_class.clone(),
        ];
        parts.extend(row.record.aliases.clone());
        parts.extend(row.record.colors.clone());
        parts.extend(row.record.assay_modes.clone());
        Self::normalize_reporter_token(&parts.join(" "))
    }

    fn normalize_reporter_constraints(
        mut constraints: ReporterConstraints,
    ) -> Result<ReporterConstraints, EngineError> {
        constraints.intended_assay = constraints
            .intended_assay
            .map(|value| Self::normalize_reporter_token(&value))
            .filter(|value| !value.is_empty());
        constraints.chassis = constraints
            .chassis
            .map(|value| Self::normalize_reporter_token(&value))
            .filter(|value| !value.is_empty());
        constraints.desired_color = constraints
            .desired_color
            .map(|value| Self::normalize_reporter_token(&value))
            .filter(|value| !value.is_empty());
        constraints.fusion_mode = constraints
            .fusion_mode
            .map(|value| Self::normalize_reporter_token(&value))
            .filter(|value| !value.is_empty());
        constraints.allowed_reporter_classes = constraints
            .allowed_reporter_classes
            .into_iter()
            .map(|value| Self::normalize_reporter_token(&value))
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let mut motifs = Vec::new();
        for motif in constraints.forbidden_motifs {
            let normalized = Self::normalize_iupac_text(&motif)?;
            if !normalized.is_empty() {
                motifs.push(normalized);
            }
        }
        constraints.forbidden_motifs = motifs;
        Ok(constraints)
    }

    fn reporter_hard_rejection_reasons(
        entry: &ReporterAnnotatedRecord,
        constraints: &ReporterConstraints,
    ) -> Result<Vec<String>, EngineError> {
        let mut reasons = Vec::new();
        let record = &entry.record;
        if !constraints.allowed_reporter_classes.is_empty()
            && !constraints
                .allowed_reporter_classes
                .contains(&Self::normalize_reporter_token(&record.reporter_class))
        {
            reasons.push(format!("class_not_allowed({})", record.reporter_class));
        }
        if let Some(max_len) = constraints.max_coding_length_bp {
            if entry.annotation.length_bp > max_len {
                reasons.push(format!(
                    "sequence_too_long({}>{})",
                    entry.annotation.length_bp, max_len
                ));
            }
        }
        if let Some(color) = constraints.desired_color.as_deref() {
            if !record
                .colors
                .iter()
                .any(|candidate| Self::normalize_reporter_token(candidate) == color)
                && Self::normalize_reporter_token(&record.spectral.color) != color
            {
                reasons.push(format!("color_mismatch({})", record.colors.join(",")));
            }
        }
        if constraints.substrate_allowed == Some(false) && record.substrate_required {
            reasons.push("substrate_required_but_not_allowed".to_string());
        }
        if constraints.live_assay == Some(true)
            && !record.assay_modes.iter().any(|mode| {
                matches!(
                    Self::normalize_reporter_token(mode).as_str(),
                    "live_cell" | "real_time"
                )
            })
        {
            reasons.push("not_marked_for_live_assay".to_string());
        }
        if !constraints.available_excitation_nm.is_empty()
            && !Self::reporter_has_channel_match(
                record.spectral.excitation_nm,
                &constraints.available_excitation_nm,
                30,
            )
        {
            reasons.push("excitation_channel_mismatch".to_string());
        }
        if !constraints.available_emission_nm.is_empty()
            && !Self::reporter_has_channel_match(
                record.spectral.emission_nm,
                &constraints.available_emission_nm,
                40,
            )
        {
            reasons.push("emission_channel_mismatch".to_string());
        }
        if let Some(fusion) = constraints.fusion_mode.as_deref() {
            let wants_fusion = matches!(fusion, "fusion" | "n_terminal" | "c_terminal");
            if wants_fusion
                && !record.fusion_compatibility.iter().any(|mode| {
                    matches!(
                        Self::normalize_reporter_token(mode).as_str(),
                        "fusion" | "n_terminal" | "c_terminal"
                    )
                })
            {
                reasons.push("not_marked_fusion_compatible".to_string());
            }
        }
        for motif in &constraints.forbidden_motifs {
            if Self::contains_motif_any_strand(record.sequence.as_bytes(), motif)? {
                reasons.push(format!("forbidden_motif_present({motif})"));
            }
        }
        Ok(reasons)
    }

    fn score_reporter_candidate(
        entry: &ReporterAnnotatedRecord,
        constraints: &ReporterConstraints,
    ) -> ReporterRecommendation {
        let record = &entry.record;
        let weights = &constraints.preference_weights;
        let mut components = BTreeMap::new();
        let confidence = match record.characterization_confidence.as_str() {
            "high" => 25.0,
            "moderate" => 15.0,
            "low" => 5.0,
            _ => 0.0,
        } * weights.characterization_confidence;
        components.insert("characterization_confidence".to_string(), confidence);

        let host_match = constraints
            .chassis
            .as_deref()
            .map(|chassis| {
                record.compatible_hosts.iter().any(|host| {
                    let host = Self::normalize_reporter_token(host);
                    host == chassis || host == "broad"
                })
            })
            .unwrap_or(false);
        let host_score = if host_match { 20.0 } else { 0.0 } * weights.host_match;
        components.insert("host_match".to_string(), host_score);

        let assay_match = constraints
            .intended_assay
            .as_deref()
            .map(|assay| {
                record
                    .assay_modes
                    .iter()
                    .any(|mode| Self::normalize_reporter_token(mode) == assay)
            })
            .unwrap_or(false);
        let assay_score = if assay_match { 20.0 } else { 0.0 } * weights.assay_match;
        components.insert("assay_match".to_string(), assay_score);

        let spectral_match = Self::reporter_has_channel_match(
            record.spectral.excitation_nm,
            &constraints.available_excitation_nm,
            30,
        ) || Self::reporter_has_channel_match(
            record.spectral.emission_nm,
            &constraints.available_emission_nm,
            40,
        );
        let spectral_score = if spectral_match { 15.0 } else { 0.0 } * weights.spectral_match;
        components.insert("spectral_match".to_string(), spectral_score);

        let brightness = record.spectral.brightness.unwrap_or(0.0).min(200.0) / 4.0;
        components.insert("brightness".to_string(), brightness * weights.brightness);

        let short_sequence = (10.0 - entry.annotation.length_bp as f64 / 300.0).max(0.0);
        components.insert(
            "short_sequence".to_string(),
            short_sequence * weights.short_sequence,
        );

        let complete_cds = if entry.annotation.likely_complete_cds {
            8.0
        } else {
            -15.0
        } * weights.complete_cds;
        components.insert("complete_cds".to_string(), complete_cds);

        let score = components.values().sum::<f64>();
        let mut rationale = Vec::new();
        if host_match {
            rationale.push("host/chassis tag matched".to_string());
        }
        if assay_match {
            rationale.push("assay mode matched".to_string());
        }
        if spectral_match {
            rationale.push("available optical channel matched".to_string());
        }
        if record.substrate_required {
            rationale.push("requires explicit substrate workflow".to_string());
        }
        if entry.annotation.likely_complete_cds {
            rationale.push("sequence passes simple complete-CDS sanity check".to_string());
        }
        if rationale.is_empty() {
            rationale.push("ranked by catalog confidence and sequence practicality".to_string());
        }

        ReporterRecommendation {
            rank: 0,
            reporter_id: record.id.clone(),
            name: record.name.clone(),
            score,
            score_components: components,
            rationale,
            warnings: entry.warnings.clone(),
            record: entry.clone(),
        }
    }

    fn reporter_has_channel_match(value: Option<u16>, channels: &[u16], tolerance_nm: u16) -> bool {
        if channels.is_empty() {
            return false;
        }
        let Some(value) = value else {
            return false;
        };
        channels
            .iter()
            .any(|channel| value.abs_diff(*channel) <= tolerance_nm)
    }

    fn normalize_reporter_token(raw: &str) -> String {
        raw.trim()
            .to_ascii_lowercase()
            .chars()
            .map(|ch| if ch.is_ascii_alphanumeric() { ch } else { '_' })
            .collect::<String>()
            .split('_')
            .filter(|part| !part.is_empty())
            .collect::<Vec<_>>()
            .join("_")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn test_record(id: &str, sequence: &str) -> ReporterRecord {
        let sequence = GentleEngine::normalize_reporter_sequence(sequence);
        ReporterRecord {
            id: id.to_string(),
            name: id.to_string(),
            reporter_class: "fluorescent_protein".to_string(),
            sequence_sha1: format!("{:x}", Sha1::digest(sequence.as_bytes())),
            sequence,
            source_refs: vec![ReporterSourceRef {
                source_id: "test".to_string(),
                accession: id.to_string(),
                url: "https://example.invalid/reporter".to_string(),
                retrieved_at: "2026-05-16".to_string(),
                license_status: "open_test".to_string(),
                license_note: "synthetic test record".to_string(),
            }],
            license_status: "open_test".to_string(),
            provenance_note: "synthetic test record".to_string(),
            colors: vec!["green".to_string()],
            assay_modes: vec!["promoter_activity".to_string(), "live_cell".to_string()],
            compatible_hosts: vec!["e_coli".to_string()],
            fusion_compatibility: vec!["standalone".to_string()],
            characterization_confidence: "high".to_string(),
            safety_scope: "benign_reporter_only".to_string(),
            spectral: ReporterSpectralProfile {
                color: "green".to_string(),
                excitation_nm: Some(501),
                emission_nm: Some(511),
                brightness: Some(50.0),
                ..ReporterSpectralProfile::default()
            },
            ..ReporterRecord::default()
        }
    }

    fn write_catalog(records: Vec<ReporterRecord>) -> (tempfile::TempDir, String) {
        let td = tempdir().expect("tempdir");
        let path = td.path().join("reporters.json");
        let catalog = ReporterCatalog {
            schema: REPORTER_CATALOG_SCHEMA.to_string(),
            curated_at: "2026-05-16".to_string(),
            records,
            ..ReporterCatalog::default()
        };
        std::fs::write(&path, serde_json::to_string_pretty(&catalog).unwrap()).unwrap();
        (td, path.to_string_lossy().to_string())
    }

    #[test]
    fn reporter_catalog_quarantines_unclear_license_and_duplicates() {
        let first = test_record("r1", "ATGAAATAA");
        let mut duplicate = test_record("r1", "ATGCCCTAA");
        duplicate.license_status = "unclear".to_string();
        duplicate.source_refs[0].license_status = "unclear".to_string();
        let (_td, path) = write_catalog(vec![first.clone(), duplicate]);
        let engine = GentleEngine::new();
        let report = engine
            .list_reporter_catalog(Some(&path), None, None)
            .expect("list catalog");
        assert_eq!(report.active_record_count, 1);
        assert_eq!(report.quarantined_record_count, 1);
        assert_eq!(report.records[0].record.id, first.id);
        assert!(report.quarantined_records[0]
            .reasons
            .contains(&"unclear_or_unaccepted_license_status".to_string()));
    }

    #[test]
    fn reporter_annotation_reports_checksum_and_orf_sanity() {
        let mut record = test_record("r1", "ATGAAATAA");
        record.sequence_sha1 = "bad".to_string();
        let annotation =
            GentleEngine::annotate_reporter_record(&record, &[]).expect("annotate reporter");
        assert_eq!(annotation.length_bp, 9);
        assert!(annotation.starts_with_atg);
        assert!(annotation.ends_with_stop);
        assert!(annotation.multiple_of_three);
        assert!(annotation.likely_complete_cds);
        assert!(!annotation.checksum_ok);
    }

    #[test]
    fn reporter_catalog_quarantines_invalid_sequence_characters() {
        let record = test_record("r1", "ATGXXXTAA");
        let (_td, path) = write_catalog(vec![record]);
        let engine = GentleEngine::new();
        let report = engine
            .list_reporter_catalog(Some(&path), None, None)
            .expect("list catalog");
        assert_eq!(report.active_record_count, 0);
        assert_eq!(report.quarantined_record_count, 1);
        assert!(report.quarantined_records[0]
            .reasons
            .contains(&"invalid_sequence_characters".to_string()));
    }

    #[test]
    fn reporter_recommendation_applies_hard_constraints_and_scores_soft_matches() {
        let mut green = test_record("green", "ATGAAATAA");
        green.colors = vec!["green".to_string()];
        green.spectral.brightness = Some(40.0);
        let mut red = test_record("red", "ATGCCCTAA");
        red.colors = vec!["red".to_string()];
        red.spectral.color = "red".to_string();
        red.spectral.excitation_nm = Some(584);
        red.spectral.emission_nm = Some(607);
        red.spectral.brightness = Some(25.0);
        let (_td, path) = write_catalog(vec![green, red]);
        let constraints = ReporterConstraints {
            desired_color: Some("red".to_string()),
            intended_assay: Some("promoter_activity".to_string()),
            chassis: Some("e_coli".to_string()),
            ..ReporterConstraints::default()
        };
        let engine = GentleEngine::new();
        let result = engine
            .recommend_reporters(Some(&path), constraints, Some(5))
            .expect("recommend");
        assert_eq!(result.recommended_candidate_count, 1);
        assert_eq!(result.recommendations[0].reporter_id, "red");
        assert_eq!(result.rejected_candidate_count, 1);
        assert!(result.rejected_candidates[0]
            .reasons
            .iter()
            .any(|reason| reason.starts_with("color_mismatch")));
    }

    #[test]
    fn reporter_corpus_jsonl_export_is_stable_and_includes_provenance() {
        let (_catalog_td, catalog_path) = write_catalog(vec![test_record("r1", "ATGAAATAA")]);
        let out_td = tempdir().expect("out tempdir");
        let output = out_td.path().join("reporters.jsonl");
        let engine = GentleEngine::new();
        let export = engine
            .export_reporter_corpus(
                Some(&catalog_path),
                &output.to_string_lossy(),
                ReporterCorpusExportFormat::Jsonl,
            )
            .expect("export corpus");
        assert_eq!(export.record_count, 1);
        let text = std::fs::read_to_string(&output).expect("read jsonl");
        assert_eq!(text.lines().count(), 1);
        assert!(text.contains("\"source_refs\""));
        assert!(text.contains("\"license_status\""));
    }
}
