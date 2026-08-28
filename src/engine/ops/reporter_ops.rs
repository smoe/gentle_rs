//! Reporter catalog validation, recommendation, and corpus export operations.

use super::*;

const DEFAULT_REPORTER_CATALOG_PATH: &str = "assets/reporter_catalog.json";
const DEFAULT_REPORTER_RECOMMENDATION_LIMIT: usize = 5;
const REPORTER_CONSTRUCT_MACRO_TEMPLATE_ID: &str = "allele_paired_promoter_luciferase_reporter";
const DEFAULT_REPORTER_BACKBONE_SEQ_ID: &str = "gentle_mammalian_luciferase_backbone_v1";
const DEFAULT_REPORTER_BACKBONE_LOAD_PATH: &str =
    "data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb";
const DEFAULT_REPORTER_CONSTRUCT_OVERLAP_BP: usize = 20;

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
        let biological_intent = Self::reporter_biological_intent(&normalized_constraints);
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
            biological_intent,
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

    pub(crate) fn plan_reporter_construct_handoff(
        &self,
        candidate_set_path: &str,
        candidate_id: Option<&str>,
        catalog_path: Option<&str>,
        reporter_constraints: Option<ReporterConstraints>,
        reporter_backbone_seq_id: Option<&str>,
        reporter_backbone_load_path: Option<&str>,
        reporter_backbone_catalog_id: Option<&str>,
        helper_catalog_path: Option<&str>,
        reference_fragment_seq_id: Option<&str>,
        alternate_fragment_seq_id: Option<&str>,
        output_prefix: Option<&str>,
    ) -> Result<ReporterConstructHandoffPlan, EngineError> {
        let candidate_set = Self::read_promoter_reporter_candidate_set(candidate_set_path)?;
        if candidate_set.schema != PROMOTER_REPORTER_CANDIDATES_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Promoter-reporter candidate set '{}' has schema '{}', expected '{}'",
                    candidate_set_path, candidate_set.schema, PROMOTER_REPORTER_CANDIDATES_SCHEMA
                ),
                cause_chain: vec![],
            });
        }
        let selected_candidate_id = candidate_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or(&candidate_set.recommended_candidate_id);
        let candidate = candidate_set
            .candidates
            .iter()
            .find(|row| row.candidate_id == selected_candidate_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Candidate '{}' was not found in promoter-reporter candidate set '{}'",
                    selected_candidate_id, candidate_set_path
                ),
                cause_chain: vec![],
            })?;

        let mut constraints =
            Self::normalize_reporter_constraints(reporter_constraints.unwrap_or_default())?;
        let luciferase_route_requested = if constraints.allowed_reporter_classes.is_empty() {
            constraints
                .allowed_reporter_classes
                .push("luciferase".to_string());
            true
        } else {
            constraints
                .allowed_reporter_classes
                .iter()
                .any(|class| Self::normalize_reporter_token(class) == "luciferase")
        };
        if constraints.intended_assay.is_none() {
            constraints.intended_assay = Some("promoter_activity".to_string());
        }
        let biological_intent = "allele_paired_promoter_luciferase_reporter_handoff".to_string();
        let reporter_biological_intent = Self::reporter_biological_intent(&constraints);

        let reporter_recommendation = if luciferase_route_requested {
            self.recommend_reporters(catalog_path, constraints.clone(), Some(1))?
        } else {
            ReporterRecommendationResult {
                schema: REPORTER_RECOMMENDATION_SCHEMA.to_string(),
                generated_at_unix_ms: Self::now_unix_ms(),
                biological_intent: reporter_biological_intent,
                catalog_path: Self::effective_reporter_catalog_label(catalog_path),
                constraints: constraints.clone(),
                warnings: vec![
                    "no_compatible_macro_route: V1 supports only luciferase reporters through allele_paired_promoter_luciferase_reporter"
                        .to_string(),
                ],
                ..ReporterRecommendationResult::default()
            }
        };
        let selected_reporter = reporter_recommendation.recommendations.first().map(|row| {
            ReporterConstructSelectedReporter {
                reporter_id: row.reporter_id.clone(),
                name: row.name.clone(),
                reporter_class: row.record.record.reporter_class.clone(),
                score: row.score,
                substrate_required: row.record.record.substrate_required,
                rationale: row.rationale.clone(),
                warnings: row.warnings.clone(),
            }
        });

        let extract_fragment_seq_id = candidate.candidate_id.clone();
        let reference_seq_id = reference_fragment_seq_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(ToOwned::to_owned)
            .unwrap_or_else(|| format!("{}_reference", extract_fragment_seq_id));
        let alternate_seq_id = alternate_fragment_seq_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(ToOwned::to_owned)
            .unwrap_or_else(|| format!("{}_alternate", extract_fragment_seq_id));
        let backbone_seq_id = reporter_backbone_seq_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or(DEFAULT_REPORTER_BACKBONE_SEQ_ID)
            .to_string();
        let backbone_load_path = reporter_backbone_load_path
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(ToOwned::to_owned)
            .or_else(|| Some(DEFAULT_REPORTER_BACKBONE_LOAD_PATH.to_string()));
        let output_prefix = output_prefix
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(ToOwned::to_owned)
            .unwrap_or_else(|| {
                format!(
                    "{}_reporter_pair",
                    Self::normalize_id_token(&candidate_set.seq_id)
                )
            });

        let selected_fragment = ReporterConstructSelectedFragment {
            candidate_id: candidate.candidate_id.clone(),
            source_seq_id: candidate_set.seq_id.clone(),
            variant_label: candidate_set.variant_label.clone(),
            gene_label: candidate.gene_label.clone(),
            transcript_id: candidate.transcript_id.clone(),
            transcript_label: candidate.transcript_label.clone(),
            start_0based: candidate.start_0based,
            end_0based_exclusive: candidate.end_0based_exclusive,
            length_bp: candidate.length_bp,
            extract_fragment_seq_id: extract_fragment_seq_id.clone(),
            reference_fragment_seq_id: reference_seq_id.clone(),
            alternate_fragment_seq_id: alternate_seq_id.clone(),
            rationale: candidate.rationale.clone(),
        };

        let backbone = self.resolve_reporter_backbone(
            &backbone_seq_id,
            backbone_load_path.as_deref(),
            reporter_backbone_seq_id.is_some(),
            reporter_backbone_catalog_id,
            helper_catalog_path,
        )?;
        let mut port_bindings = vec![
            self.reporter_sequence_port_binding(
                "reference_fragment_seq_id",
                &reference_seq_id,
                reference_fragment_seq_id.is_some(),
                true,
                "Reference-allele promoter fragment for the macro template.",
            ),
            self.reporter_sequence_port_binding(
                "alternate_fragment_seq_id",
                &alternate_seq_id,
                alternate_fragment_seq_id.is_some(),
                true,
                "Alternate-allele promoter fragment for the macro template.",
            ),
            self.reporter_backbone_port_binding(&backbone, reporter_backbone_seq_id.is_some()),
            ReporterConstructPortBinding {
                port_id: "overlap_bp".to_string(),
                value: Some(DEFAULT_REPORTER_CONSTRUCT_OVERLAP_BP.to_string()),
                status: PortBindingStatus::Ready,
                required: true,
                note: "Stable reserved overlap parameter for the existing macro template."
                    .to_string(),
            },
        ];
        if !luciferase_route_requested || selected_reporter.is_none() {
            port_bindings.clear();
        }

        let mut warnings = reporter_recommendation.warnings.clone();
        if !luciferase_route_requested {
            warnings.push(
                "no_compatible_macro_route: V1 supports only luciferase reporters through allele_paired_promoter_luciferase_reporter"
                    .to_string(),
            );
        } else if selected_reporter.is_none() {
            warnings.push(
                "No luciferase reporter candidate satisfied the handoff constraints".to_string(),
            );
        }

        let exact_vector_verified = backbone
            .validation
            .as_ref()
            .map(|report| report.status == ReporterVectorValidationStatus::Verified)
            .unwrap_or(true);
        if let Some(validation) = backbone.validation.as_ref() {
            warnings.extend(validation.warnings.iter().cloned());
            if validation.status != ReporterVectorValidationStatus::Verified {
                warnings.push(format!(
                    "Exact reporter vector '{}' is {}; construct commands remain disabled until the catalog identity is verified",
                    validation.helper_catalog_id,
                    match validation.status {
                        ReporterVectorValidationStatus::Verified => "verified",
                        ReporterVectorValidationStatus::Rejected => "rejected",
                        ReporterVectorValidationStatus::Unavailable => "unavailable",
                    }
                ));
            }
        }

        let commands =
            if luciferase_route_requested && selected_reporter.is_some() && exact_vector_verified {
                Self::reporter_construct_handoff_commands(
                    &candidate_set,
                    &candidate,
                    &extract_fragment_seq_id,
                    &reference_seq_id,
                    &alternate_seq_id,
                    &backbone,
                    &output_prefix,
                )
            } else {
                vec![]
            };
        let status = if !luciferase_route_requested {
            "no_compatible_macro_route"
        } else if selected_reporter.is_none() {
            "no_compatible_reporter"
        } else if backbone.status == ReporterBackboneResolutionStatus::ExactIdentityRejected {
            "backbone_validation_rejected"
        } else if backbone.status == ReporterBackboneResolutionStatus::ExactIdentityUnavailable {
            "backbone_validation_unavailable"
        } else if port_bindings
            .iter()
            .filter(|binding| binding.required)
            .all(|binding| binding.status == PortBindingStatus::Ready)
        {
            "ready"
        } else {
            "needs_inputs"
        }
        .to_string();

        Ok(ReporterConstructHandoffPlan {
            schema: REPORTER_CONSTRUCT_HANDOFF_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            status,
            biological_intent,
            provenance: ReporterConstructHandoffProvenance {
                candidate_set_path: candidate_set_path.to_string(),
                candidate_set_schema: candidate_set.schema.clone(),
                candidate_set_generated_at_unix_ms: candidate_set.generated_at_unix_ms,
                candidate_set_op_id: candidate_set.op_id.clone(),
                candidate_set_run_id: candidate_set.run_id.clone(),
                reporter_catalog_path: reporter_recommendation.catalog_path.clone(),
                macro_template_id: REPORTER_CONSTRUCT_MACRO_TEMPLATE_ID.to_string(),
            },
            selected_fragment,
            selected_reporter,
            reporter_recommendation,
            backbone,
            port_bindings,
            commands,
            warnings,
        })
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
        if record.sequence_sha1.trim().is_empty()
            && record
                .sequence_sha256
                .as_deref()
                .unwrap_or_default()
                .trim()
                .is_empty()
        {
            reasons.push("missing_sequence_checksum".to_string());
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
        let multiple_of_three = sequence.len().is_multiple_of(3);
        let checksum_sha256 = crate::digest_utils::sha256_hex_str(&record.sequence);
        let checksum_ok = if let Some(expected) = record
            .sequence_sha256
            .as_deref()
            .and_then(Self::normalize_reporter_sha256_checksum)
        {
            expected == checksum_sha256
        } else {
            Self::looks_like_legacy_sha1_checksum(&record.sequence_sha1)
        };
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
            checksum_ok,
            forbidden_motif_hits,
        })
    }

    fn normalize_reporter_sha256_checksum(raw: &str) -> Option<String> {
        let trimmed = raw.trim();
        let value = trimmed.strip_prefix("sha256:").unwrap_or(trimmed);
        if value.len() == 64 && value.as_bytes().iter().all(u8::is_ascii_hexdigit) {
            Some(value.to_ascii_lowercase())
        } else {
            None
        }
    }

    fn looks_like_legacy_sha1_checksum(raw: &str) -> bool {
        let trimmed = raw.trim();
        trimmed.len() == 40 && trimmed.as_bytes().iter().all(u8::is_ascii_hexdigit)
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

    fn reporter_biological_intent(constraints: &ReporterConstraints) -> String {
        let wants_luciferase = constraints
            .allowed_reporter_classes
            .iter()
            .any(|class| class == "luciferase");
        let assay = constraints.intended_assay.as_deref().unwrap_or_default();
        let promoter_like = assay.contains("promoter") || assay.contains("regulatory");
        if promoter_like && wants_luciferase {
            return "promoter_luciferase_reporter_assay".to_string();
        }
        if promoter_like {
            return "promoter_reporter_assay".to_string();
        }
        if let Some(fusion) = constraints.fusion_mode.as_deref()
            && matches!(fusion, "fusion" | "n_terminal" | "c_terminal")
        {
            return "fusion_reporter_assay".to_string();
        }
        if constraints.live_assay == Some(true) {
            return "live_reporter_assay".to_string();
        }
        if constraints.live_assay == Some(false) {
            return "endpoint_reporter_assay".to_string();
        }
        if wants_luciferase {
            return "luciferase_reporter_selection".to_string();
        }
        "reporter_selection".to_string()
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
        if let Some(max_len) = constraints.max_coding_length_bp
            && entry.annotation.length_bp > max_len
        {
            reasons.push(format!(
                "sequence_too_long({}>{})",
                entry.annotation.length_bp, max_len
            ));
        }
        if let Some(color) = constraints.desired_color.as_deref()
            && !record
                .colors
                .iter()
                .any(|candidate| Self::normalize_reporter_token(candidate) == color)
            && Self::normalize_reporter_token(&record.spectral.color) != color
        {
            reasons.push(format!("color_mismatch({})", record.colors.join(",")));
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

    fn read_promoter_reporter_candidate_set(
        path: &str,
    ) -> Result<PromoterReporterCandidateSet, EngineError> {
        let text = std::fs::read_to_string(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read promoter-reporter candidate set '{path}': {e}"),
            cause_chain: vec![],
        })?;
        serde_json::from_str(&text).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not parse promoter-reporter candidate set '{path}': {e}"),
            cause_chain: vec![],
        })
    }

    fn resolve_reporter_backbone(
        &self,
        seq_id: &str,
        load_path: Option<&str>,
        explicit_seq_id: bool,
        helper_catalog_id: Option<&str>,
        helper_catalog_path: Option<&str>,
    ) -> Result<ReporterBackboneResolution, EngineError> {
        let exact_catalog_id = helper_catalog_id
            .map(str::trim)
            .filter(|value| !value.is_empty());
        if let Some(helper_catalog_id) = exact_catalog_id {
            let (catalog, catalog_origin) = Self::open_helper_genome_catalog(helper_catalog_path)?;
            let expectation = catalog
                .helper_vector_sequence_expectation(helper_catalog_id)
                .map_err(|message| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Could not resolve exact reporter vector '{}' from helper catalog '{}': {}",
                        helper_catalog_id, catalog_origin, message
                    ),
                    cause_chain: vec![],
                })?
                .ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Helper catalog entry '{}' has no sequence_expectation and cannot gate exact reporter-vector identity",
                        helper_catalog_id
                    ),
                    cause_chain: vec![],
                })?;

            let loaded_from_path = if self.state.sequences.contains_key(seq_id) {
                None
            } else if let Some(path) = load_path {
                match crate::dna_sequence::load_from_file(path) {
                    Ok(mut dna) => {
                        Self::prepare_sequence(&mut dna);
                        Some(Ok(dna))
                    }
                    Err(error) => Some(Err(format!(
                        "Could not load reporter backbone path '{}' for exact validation: {}",
                        path, error
                    ))),
                }
            } else {
                None
            };
            let validation = if let Some(dna) = self.state.sequences.get(seq_id) {
                Self::validate_reporter_vector_sequence(
                    seq_id,
                    helper_catalog_id,
                    &catalog_origin,
                    &expectation,
                    dna,
                )
            } else if let Some(Ok(dna)) = loaded_from_path.as_ref() {
                Self::validate_reporter_vector_sequence(
                    seq_id,
                    helper_catalog_id,
                    &catalog_origin,
                    &expectation,
                    dna,
                )
            } else {
                let warning = loaded_from_path.and_then(Result::err).unwrap_or_else(|| {
                    format!(
                        "No loaded sequence or readable backbone path was available for '{}'",
                        seq_id
                    )
                });
                Self::unavailable_reporter_vector_validation(
                    seq_id,
                    helper_catalog_id,
                    &catalog_origin,
                    &expectation,
                    warning,
                )
            };
            let status = match validation.status {
                ReporterVectorValidationStatus::Verified
                    if self.state.sequences.contains_key(seq_id) =>
                {
                    ReporterBackboneResolutionStatus::ResolvedInState
                }
                ReporterVectorValidationStatus::Verified => {
                    ReporterBackboneResolutionStatus::RequiresManualLoad
                }
                ReporterVectorValidationStatus::Rejected => {
                    ReporterBackboneResolutionStatus::ExactIdentityRejected
                }
                ReporterVectorValidationStatus::Unavailable => {
                    ReporterBackboneResolutionStatus::ExactIdentityUnavailable
                }
            };
            let note = match status {
                ReporterBackboneResolutionStatus::ResolvedInState => {
                    "Exact reporter-vector identity is verified for the sequence in project state."
                }
                ReporterBackboneResolutionStatus::RequiresManualLoad => {
                    "Exact reporter-vector identity is verified at the supplied path; load it before materialization."
                }
                ReporterBackboneResolutionStatus::ExactIdentityRejected => {
                    "The supplied sequence failed exact reporter-vector identity validation; construct commands are disabled."
                }
                ReporterBackboneResolutionStatus::ExactIdentityUnavailable => {
                    "Exact reporter-vector identity could not be evaluated; construct commands are disabled."
                }
                ReporterBackboneResolutionStatus::UnresolvedSeqIdProvided => unreachable!(),
            };
            return Ok(ReporterBackboneResolution {
                seq_id: seq_id.to_string(),
                load_path: load_path.map(ToOwned::to_owned),
                status,
                note: note.to_string(),
                helper_catalog_id: Some(helper_catalog_id.to_string()),
                validation: Some(validation),
            });
        }

        if self.state.sequences.contains_key(seq_id) {
            return Ok(ReporterBackboneResolution {
                seq_id: seq_id.to_string(),
                load_path: load_path.map(ToOwned::to_owned),
                status: ReporterBackboneResolutionStatus::ResolvedInState,
                note: "Reporter backbone sequence is already loaded in project state.".to_string(),
                helper_catalog_id: None,
                validation: None,
            });
        }
        if load_path.is_some() {
            return Ok(ReporterBackboneResolution {
                seq_id: seq_id.to_string(),
                load_path: load_path.map(ToOwned::to_owned),
                status: ReporterBackboneResolutionStatus::RequiresManualLoad,
                note: "Reporter backbone must be loaded before running the macro template."
                    .to_string(),
                helper_catalog_id: None,
                validation: None,
            });
        }
        Ok(ReporterBackboneResolution {
            seq_id: seq_id.to_string(),
            load_path: None,
            status: ReporterBackboneResolutionStatus::UnresolvedSeqIdProvided,
            note: if explicit_seq_id {
                "Reporter backbone seq_id was provided but is not present in project state."
                    .to_string()
            } else {
                "Reporter backbone seq_id is not present and no load path was provided.".to_string()
            },
            helper_catalog_id: None,
            validation: None,
        })
    }

    pub(crate) fn validate_reporter_vector_sequence(
        seq_id: &str,
        helper_catalog_id: &str,
        helper_catalog_path: &str,
        expectation: &HelperVectorSequenceExpectation,
        dna: &DNAsequence,
    ) -> ReporterVectorValidationReport {
        let observed_accession = dna.accession().map(str::to_string);
        let observed_version = dna.version().map(str::to_string);
        let observed_identity = observed_version
            .as_deref()
            .or(observed_accession.as_deref())
            .unwrap_or("missing");
        let observed_topology = if dna.is_circular() {
            "circular"
        } else {
            "linear"
        };
        let mut checks = vec![
            ReporterVectorValidationCheck {
                check_id: "accession_version".to_string(),
                required: true,
                passed: observed_identity.eq_ignore_ascii_case(&expectation.accession_version),
                expected: expectation.accession_version.clone(),
                observed: observed_identity.to_string(),
                detail: "The versioned GenBank identity must match exactly.".to_string(),
            },
            ReporterVectorValidationCheck {
                check_id: "length_bp".to_string(),
                required: true,
                passed: dna.len() == expectation.expected_length_bp,
                expected: expectation.expected_length_bp.to_string(),
                observed: dna.len().to_string(),
                detail: "Sequence length is catalog-owned exact identity evidence.".to_string(),
            },
            ReporterVectorValidationCheck {
                check_id: "topology".to_string(),
                required: true,
                passed: observed_topology.eq_ignore_ascii_case(&expectation.expected_topology),
                expected: expectation.expected_topology.clone(),
                observed: observed_topology.to_string(),
                detail: "Reporter-vector topology must match the catalog expectation.".to_string(),
            },
        ];

        let mut observed_luc2_start_1based = None;
        let mut observed_multiple_cloning_region = None;
        for required_feature in &expectation.required_features {
            let matched = Self::find_reporter_vector_feature(dna, required_feature);
            let expected_interval = Self::format_expected_feature_interval(required_feature);
            let (passed, observed, detail) = if let Some((start_1based, end_1based, feature_text)) =
                matched
            {
                let coordinates_match = required_feature
                    .expected_start_1based
                    .is_none_or(|expected| expected == start_1based)
                    && required_feature
                        .expected_end_1based
                        .is_none_or(|expected| expected == end_1based);
                if required_feature.id.eq_ignore_ascii_case("luc2") {
                    observed_luc2_start_1based = Some(start_1based);
                }
                if required_feature
                    .id
                    .eq_ignore_ascii_case("multiple_cloning_region")
                {
                    observed_multiple_cloning_region =
                        Some(format!("{}..{}", start_1based, end_1based));
                }
                (
                    coordinates_match,
                    format!("{}..{} ({})", start_1based, end_1based, feature_text),
                    if coordinates_match {
                        "Required annotation and boundary matched.".to_string()
                    } else {
                        "Required annotation matched, but its boundary differs.".to_string()
                    },
                )
            } else {
                (
                    false,
                    "missing".to_string(),
                    "No annotated feature matched the required kind/qualifier terms.".to_string(),
                )
            };
            checks.push(ReporterVectorValidationCheck {
                check_id: format!("feature:{}", required_feature.id),
                required: true,
                passed,
                expected: expected_interval,
                observed,
                detail,
            });
        }

        let mut site_counts: BTreeMap<String, usize> = BTreeMap::new();
        for site in dna
            .restriction_enzyme_sites()
            .iter()
            .filter(|site| site.forward_strand)
        {
            *site_counts.entry(site.enzyme.name.clone()).or_default() += 1;
        }
        let unique_restriction_sites = site_counts
            .into_iter()
            .filter_map(|(name, count)| (count == 1).then_some(name))
            .collect::<Vec<_>>();
        checks.push(ReporterVectorValidationCheck {
            check_id: "derived_unique_restriction_sites".to_string(),
            required: false,
            passed: true,
            expected: "derived from the imported sequence and active enzyme catalog".to_string(),
            observed: unique_restriction_sites.join(","),
            detail:
                "This inventory is reported for cloning review and is not an identity assertion."
                    .to_string(),
        });

        let status = if checks
            .iter()
            .filter(|check| check.required)
            .all(|check| check.passed)
        {
            ReporterVectorValidationStatus::Verified
        } else {
            ReporterVectorValidationStatus::Rejected
        };
        let expected_luc2_start_1based = expectation
            .required_features
            .iter()
            .find(|feature| feature.id.eq_ignore_ascii_case("luc2"))
            .and_then(|feature| feature.expected_start_1based);
        let expected_multiple_cloning_region = expectation
            .required_features
            .iter()
            .find(|feature| feature.id.eq_ignore_ascii_case("multiple_cloning_region"))
            .map(Self::format_expected_feature_interval);
        ReporterVectorValidationReport {
            schema: REPORTER_VECTOR_VALIDATION_SCHEMA.to_string(),
            helper_catalog_id: helper_catalog_id.to_string(),
            helper_catalog_path: helper_catalog_path.to_string(),
            seq_id: seq_id.to_string(),
            status,
            expected_accession_version: expectation.accession_version.clone(),
            observed_accession,
            observed_version,
            expected_length_bp: expectation.expected_length_bp,
            observed_length_bp: Some(dna.len()),
            expected_topology: expectation.expected_topology.clone(),
            observed_topology: Some(observed_topology.to_string()),
            expected_luc2_start_1based,
            observed_luc2_start_1based,
            expected_multiple_cloning_region,
            observed_multiple_cloning_region,
            unique_restriction_sites,
            restriction_site_equivalence_notes: expectation
                .restriction_site_equivalences
                .iter()
                .map(|row| {
                    format!(
                        "{} / {}: {}",
                        row.detected_enzyme, row.equivalent_enzyme, row.note
                    )
                })
                .collect(),
            checks,
            warnings: vec![
                "Exact sequence identity validation supports construct planning; it does not verify physical stock identity."
                    .to_string(),
            ],
        }
    }

    pub(crate) fn validate_known_helper_vector_accession(
        &self,
        seq_id: &str,
        accession: &str,
    ) -> Result<Option<ReporterVectorValidationReport>, EngineError> {
        let (catalog, catalog_origin) = Self::open_helper_genome_catalog(None)?;
        let Some((helper_catalog_id, expectation)) =
            catalog.helper_vector_sequence_expectation_for_accession(accession)
        else {
            return Ok(None);
        };
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Fetched helper-vector sequence '{}' disappeared before validation",
                    seq_id
                ),
                cause_chain: vec![],
            })?;
        Ok(Some(Self::validate_reporter_vector_sequence(
            seq_id,
            &helper_catalog_id,
            &catalog_origin,
            &expectation,
            dna,
        )))
    }

    fn unavailable_reporter_vector_validation(
        seq_id: &str,
        helper_catalog_id: &str,
        helper_catalog_path: &str,
        expectation: &HelperVectorSequenceExpectation,
        warning: String,
    ) -> ReporterVectorValidationReport {
        let mut checks = vec![
            ("accession_version", expectation.accession_version.clone()),
            ("length_bp", expectation.expected_length_bp.to_string()),
            ("topology", expectation.expected_topology.clone()),
        ]
        .into_iter()
        .map(|(check_id, expected)| ReporterVectorValidationCheck {
            check_id: check_id.to_string(),
            required: true,
            passed: false,
            expected,
            observed: "unavailable".to_string(),
            detail: "No sequence was available for this required check.".to_string(),
        })
        .collect::<Vec<_>>();
        checks.extend(expectation.required_features.iter().map(|feature| {
            ReporterVectorValidationCheck {
                check_id: format!("feature:{}", feature.id),
                required: true,
                passed: false,
                expected: Self::format_expected_feature_interval(feature),
                observed: "unavailable".to_string(),
                detail: "No sequence annotation was available for this required check.".to_string(),
            }
        }));
        ReporterVectorValidationReport {
            schema: REPORTER_VECTOR_VALIDATION_SCHEMA.to_string(),
            helper_catalog_id: helper_catalog_id.to_string(),
            helper_catalog_path: helper_catalog_path.to_string(),
            seq_id: seq_id.to_string(),
            status: ReporterVectorValidationStatus::Unavailable,
            expected_accession_version: expectation.accession_version.clone(),
            expected_length_bp: expectation.expected_length_bp,
            expected_topology: expectation.expected_topology.clone(),
            expected_luc2_start_1based: expectation
                .required_features
                .iter()
                .find(|feature| feature.id.eq_ignore_ascii_case("luc2"))
                .and_then(|feature| feature.expected_start_1based),
            expected_multiple_cloning_region: expectation
                .required_features
                .iter()
                .find(|feature| feature.id.eq_ignore_ascii_case("multiple_cloning_region"))
                .map(Self::format_expected_feature_interval),
            restriction_site_equivalence_notes: expectation
                .restriction_site_equivalences
                .iter()
                .map(|row| {
                    format!(
                        "{} / {}: {}",
                        row.detected_enzyme, row.equivalent_enzyme, row.note
                    )
                })
                .collect(),
            checks,
            warnings: vec![warning],
            ..ReporterVectorValidationReport::default()
        }
    }

    fn find_reporter_vector_feature(
        dna: &DNAsequence,
        expectation: &HelperVectorRequiredFeatureExpectation,
    ) -> Option<(usize, usize, String)> {
        let mut matches = dna
            .features()
            .iter()
            .filter_map(|feature| {
                let kind = feature.kind.to_string();
                let kind_matches = expectation.feature_kinds.is_empty()
                    || expectation
                        .feature_kinds
                        .iter()
                        .any(|expected| kind.eq_ignore_ascii_case(expected));
                if !kind_matches {
                    return None;
                }
                let mut text_parts = vec![kind.clone()];
                for (key, value) in &feature.qualifiers {
                    let key = key.to_string();
                    if matches!(key.as_str(), "label" | "gene" | "product" | "note") {
                        text_parts.push(match value.as_deref() {
                            Some(value) => format!("{key}={value}"),
                            None => key,
                        });
                    }
                }
                let text = text_parts.join(" ");
                let lower = text.to_ascii_lowercase();
                let terms_match = expectation.qualifier_terms.is_empty()
                    || expectation
                        .qualifier_terms
                        .iter()
                        .any(|term| lower.contains(&term.trim().to_ascii_lowercase()));
                if !terms_match {
                    return None;
                }
                let (start, end) = feature.location.find_bounds().ok()?;
                let start = usize::try_from(start).ok()?.saturating_add(1);
                let end = usize::try_from(end).ok()?;
                Some((start, end, text))
            })
            .collect::<Vec<_>>();
        matches.sort_by(|left, right| {
            (left.0, left.1, left.2.as_str()).cmp(&(right.0, right.1, right.2.as_str()))
        });
        matches.into_iter().next()
    }

    fn format_expected_feature_interval(
        expectation: &HelperVectorRequiredFeatureExpectation,
    ) -> String {
        let interval = match (
            expectation.expected_start_1based,
            expectation.expected_end_1based,
        ) {
            (Some(start), Some(end)) => format!("{}..{}", start, end),
            (Some(start), None) => format!("start={}", start),
            (None, Some(end)) => format!("end={}", end),
            (None, None) => "annotated feature".to_string(),
        };
        format!(
            "{} [{}; terms={}]",
            interval,
            expectation.feature_kinds.join("|"),
            expectation.qualifier_terms.join("|")
        )
    }

    fn reporter_sequence_port_binding(
        &self,
        port_id: &str,
        seq_id: &str,
        explicit_seq_id: bool,
        derivable: bool,
        note: &str,
    ) -> ReporterConstructPortBinding {
        let status = if self.state.sequences.contains_key(seq_id) {
            PortBindingStatus::Ready
        } else if explicit_seq_id {
            PortBindingStatus::ProvidedMissingFromState
        } else if derivable {
            PortBindingStatus::Derivable
        } else {
            PortBindingStatus::Missing
        };
        ReporterConstructPortBinding {
            port_id: port_id.to_string(),
            value: Some(seq_id.to_string()),
            status,
            required: true,
            note: note.to_string(),
        }
    }

    fn reporter_backbone_port_binding(
        &self,
        backbone: &ReporterBackboneResolution,
        explicit_seq_id: bool,
    ) -> ReporterConstructPortBinding {
        let status = match backbone.status {
            ReporterBackboneResolutionStatus::ResolvedInState => PortBindingStatus::Ready,
            ReporterBackboneResolutionStatus::UnresolvedSeqIdProvided => {
                if explicit_seq_id {
                    PortBindingStatus::ProvidedMissingFromState
                } else {
                    PortBindingStatus::Missing
                }
            }
            ReporterBackboneResolutionStatus::RequiresManualLoad => {
                if explicit_seq_id {
                    PortBindingStatus::ProvidedMissingFromState
                } else {
                    PortBindingStatus::Missing
                }
            }
            ReporterBackboneResolutionStatus::ExactIdentityUnavailable
            | ReporterBackboneResolutionStatus::ExactIdentityRejected => PortBindingStatus::Missing,
        };
        ReporterConstructPortBinding {
            port_id: "reporter_backbone_seq_id".to_string(),
            value: Some(backbone.seq_id.clone()),
            status,
            required: true,
            note: backbone.note.clone(),
        }
    }

    fn reporter_construct_handoff_commands(
        candidate_set: &PromoterReporterCandidateSet,
        candidate: &PromoterReporterFragmentCandidate,
        extract_fragment_seq_id: &str,
        reference_seq_id: &str,
        alternate_seq_id: &str,
        backbone: &ReporterBackboneResolution,
        output_prefix: &str,
    ) -> Vec<ReporterConstructHandoffCommand> {
        let mut commands = vec![ReporterConstructHandoffCommand {
            label: "Extract selected promoter fragment".to_string(),
            command_kind: "op".to_string(),
            command: format!(
                "op {}",
                serde_json::json!({
                    "ExtractRegion": {
                        "input": candidate_set.seq_id,
                        "from": candidate.start_0based,
                        "to": candidate.end_0based_exclusive,
                        "output_id": extract_fragment_seq_id,
                    }
                })
            ),
            mutating: true,
            note: "Creates the selected promoter fragment sequence before allele materialization."
                .to_string(),
        }];
        commands.push(ReporterConstructHandoffCommand {
            label: "Materialize reference allele fragment".to_string(),
            command_kind: "op".to_string(),
            command: format!(
                "op {}",
                serde_json::json!({
                    "MaterializeVariantAllele": {
                        "input": extract_fragment_seq_id,
                        "variant_label_or_id": candidate_set.variant_label,
                        "allele": "reference",
                        "output_id": reference_seq_id,
                    }
                })
            ),
            mutating: true,
            note: "Creates the reference-allele promoter insert for the macro template."
                .to_string(),
        });
        commands.push(ReporterConstructHandoffCommand {
            label: "Materialize alternate allele fragment".to_string(),
            command_kind: "op".to_string(),
            command: format!(
                "op {}",
                serde_json::json!({
                    "MaterializeVariantAllele": {
                        "input": extract_fragment_seq_id,
                        "variant_label_or_id": candidate_set.variant_label,
                        "allele": "alternate",
                        "output_id": alternate_seq_id,
                    }
                })
            ),
            mutating: true,
            note: "Creates the alternate-allele promoter insert for the macro template."
                .to_string(),
        });
        if backbone.status == ReporterBackboneResolutionStatus::RequiresManualLoad
            && let Some(load_path) = backbone.load_path.as_deref()
        {
            commands.push(ReporterConstructHandoffCommand {
                label: "Load reporter backbone".to_string(),
                command_kind: "op".to_string(),
                command: format!(
                    "op {}",
                    serde_json::json!({
                        "LoadFile": {
                            "path": load_path,
                            "as_id": backbone.seq_id,
                        }
                    })
                ),
                mutating: true,
                note: "Loads the promoterless luciferase reporter backbone into project state."
                    .to_string(),
            });
        }
        commands.push(ReporterConstructHandoffCommand {
            label: "Import reporter macro templates".to_string(),
            command_kind: "shell".to_string(),
            command: "macros template-import assets/cloning_patterns_catalog".to_string(),
            mutating: true,
            note: "Makes the built-in reporter macro template available in the project."
                .to_string(),
        });
        let macro_bindings = format!(
            "macros template-run {} --bind reference_fragment_seq_id={} --bind alternate_fragment_seq_id={} --bind reporter_backbone_seq_id={} --bind overlap_bp={} --bind output_prefix={} --transactional",
            REPORTER_CONSTRUCT_MACRO_TEMPLATE_ID,
            Self::quote_reporter_shell_token(reference_seq_id),
            Self::quote_reporter_shell_token(alternate_seq_id),
            Self::quote_reporter_shell_token(&backbone.seq_id),
            DEFAULT_REPORTER_CONSTRUCT_OVERLAP_BP,
            Self::quote_reporter_shell_token(output_prefix)
        );
        commands.push(ReporterConstructHandoffCommand {
            label: "Validate reporter macro".to_string(),
            command_kind: "shell".to_string(),
            command: format!("{macro_bindings} --validate-only"),
            mutating: false,
            note: "Checks macro-template bindings before creating construct preview sequences."
                .to_string(),
        });
        commands.push(ReporterConstructHandoffCommand {
            label: "Run reporter macro".to_string(),
            command_kind: "shell".to_string(),
            command: macro_bindings,
            mutating: true,
            note: "Creates paired reference/alternate promoter-luciferase construct previews."
                .to_string(),
        });
        commands
    }

    fn quote_reporter_shell_token(value: &str) -> String {
        if !value.is_empty()
            && value
                .chars()
                .all(|ch| ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.' | '/' | ':'))
        {
            return value.to_string();
        }
        format!("'{}'", value.replace('\'', "'\\''"))
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
            sequence_sha1: String::new(),
            sequence_sha256: Some(crate::digest_utils::sha256_hex_str(&sequence)),
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

    fn vkorc1_candidate_set_path() -> &'static str {
        "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/promoter_reporter_candidates.json"
    }

    fn synthetic_mcs_expectation() -> HelperVectorSequenceExpectation {
        HelperVectorSequenceExpectation {
            schema: crate::genomes::HELPER_VECTOR_SEQUENCE_EXPECTATION_SCHEMA.to_string(),
            provider: "GENtle tests".to_string(),
            product_name: "synthetic MCS backbone".to_string(),
            catalog_number: "SYNTH-MCS-1".to_string(),
            accession_version: "GENTLE_SYNTHETIC_MCS.1".to_string(),
            expected_length_bp: 240,
            expected_topology: "circular".to_string(),
            required_features: vec![
                HelperVectorRequiredFeatureExpectation {
                    id: "multiple_cloning_region".to_string(),
                    feature_kinds: vec!["misc_feature".to_string()],
                    qualifier_terms: vec!["multiple cloning site region".to_string()],
                    expected_start_1based: Some(1),
                    expected_end_1based: Some(70),
                },
                HelperVectorRequiredFeatureExpectation {
                    id: "luc2".to_string(),
                    feature_kinds: vec!["CDS".to_string()],
                    qualifier_terms: vec!["luciferase luc2 marker".to_string()],
                    expected_start_1based: Some(100),
                    expected_end_1based: Some(180),
                },
            ],
            restriction_site_equivalences: vec![
                crate::genomes::HelperVectorRestrictionSiteEquivalence {
                    detected_enzyme: "KpnI".to_string(),
                    equivalent_enzyme: "Acc65I".to_string(),
                    note: "Synthetic fixture uses the shared recognition site.".to_string(),
                },
                crate::genomes::HelperVectorRestrictionSiteEquivalence {
                    detected_enzyme: "SacI".to_string(),
                    equivalent_enzyme: "EcoICRI".to_string(),
                    note: "Synthetic fixture uses the shared recognition site.".to_string(),
                },
            ],
            provenance: vec![crate::genomes::HelperVectorSequenceExpectationProvenance {
                source_id: "synthetic-test-fixture".to_string(),
                source_url: "test_files/fixtures/reporter_vectors/synthetic_mcs_backbone.gb"
                    .to_string(),
                asserted_on: "2026-08-27".to_string(),
                note: "Repository-owned deterministic fixture.".to_string(),
            }],
            ..HelperVectorSequenceExpectation::default()
        }
    }

    #[test]
    fn pgl4_catalog_entry_carries_exact_versioned_identity() {
        let catalog = GenomeCatalog::from_json_file("assets/helper_genomes.json")
            .expect("load helper catalog");
        let expectation = catalog
            .helper_vector_sequence_expectation("Reporter Promega Luciferase AY738222 (online)")
            .expect("resolve pGL4 helper")
            .expect("pGL4 exact expectation");

        assert_eq!(expectation.provider, "Promega");
        assert_eq!(expectation.product_name, "pGL4.10[luc2]");
        assert_eq!(expectation.catalog_number, "E6651");
        assert_eq!(expectation.accession_version, "AY738222.1");
        assert_eq!(expectation.expected_length_bp, 4242);
        assert_eq!(expectation.expected_topology, "circular");
        assert!(
            expectation
                .required_features
                .iter()
                .any(|feature| feature.id == "luc2" && feature.expected_start_1based == Some(100))
        );
    }

    #[test]
    fn synthetic_mcs_fixture_validates_against_its_own_identity() {
        let mut dna = crate::dna_sequence::load_from_file(
            "test_files/fixtures/reporter_vectors/synthetic_mcs_backbone.gb",
        )
        .expect("load synthetic MCS fixture");
        GentleEngine::prepare_sequence(&mut dna);
        let report = GentleEngine::validate_reporter_vector_sequence(
            "synthetic_mcs",
            "Synthetic MCS fixture",
            "inline-test-expectation",
            &synthetic_mcs_expectation(),
            &dna,
        );

        assert_eq!(report.status, ReporterVectorValidationStatus::Verified);
        assert_eq!(report.observed_luc2_start_1based, Some(100));
        assert_eq!(
            report.observed_multiple_cloning_region.as_deref(),
            Some("1..70")
        );
        assert!(
            report
                .unique_restriction_sites
                .contains(&"KpnI".to_string())
        );
        assert!(
            report
                .unique_restriction_sites
                .contains(&"SacI".to_string())
        );
        assert!(
            !report
                .unique_restriction_sites
                .contains(&"Acc65I".to_string())
        );
        assert!(
            !report
                .unique_restriction_sites
                .contains(&"EcoICRI".to_string())
        );
        assert!(
            report
                .restriction_site_equivalence_notes
                .iter()
                .any(|note| note.contains("KpnI / Acc65I"))
        );
        assert!(
            report
                .restriction_site_equivalence_notes
                .iter()
                .any(|note| note.contains("SacI / EcoICRI"))
        );
    }

    #[test]
    fn synthetic_tutorial_backbone_is_rejected_as_exact_pgl4() {
        let engine = GentleEngine::new();
        let plan = engine
            .plan_reporter_construct_handoff(
                vkorc1_candidate_set_path(),
                None,
                None,
                None,
                Some("claimed_pgl4"),
                Some(DEFAULT_REPORTER_BACKBONE_LOAD_PATH),
                Some("Reporter Promega Luciferase AY738222 (online)"),
                Some("assets/helper_genomes.json"),
                None,
                None,
                None,
            )
            .expect("return a fail-closed handoff report");

        assert_eq!(plan.status, "backbone_validation_rejected");
        assert_eq!(
            plan.backbone.status,
            ReporterBackboneResolutionStatus::ExactIdentityRejected
        );
        let validation = plan.backbone.validation.expect("validation report");
        assert_eq!(validation.status, ReporterVectorValidationStatus::Rejected);
        assert_eq!(validation.expected_length_bp, 4242);
        assert_eq!(validation.observed_length_bp, Some(3301));
        assert!(plan.commands.is_empty());
    }

    #[test]
    fn exact_pgl4_request_without_sequence_is_unavailable_and_emits_no_commands() {
        let td = tempdir().expect("tempdir");
        let missing = td.path().join("not-downloaded.gb");
        let engine = GentleEngine::new();
        let plan = engine
            .plan_reporter_construct_handoff(
                vkorc1_candidate_set_path(),
                None,
                None,
                None,
                Some("pgl4_exact"),
                Some(&missing.to_string_lossy()),
                Some("Reporter Promega Luciferase AY738222 (online)"),
                Some("assets/helper_genomes.json"),
                None,
                None,
                None,
            )
            .expect("return unavailable exact-vector report");

        assert_eq!(plan.status, "backbone_validation_unavailable");
        assert_eq!(
            plan.backbone.status,
            ReporterBackboneResolutionStatus::ExactIdentityUnavailable
        );
        assert!(plan.commands.is_empty());
    }

    #[test]
    fn real_pgl4_genbank_validation_is_online_opt_in() {
        if std::env::var_os("GENTLE_TEST_ONLINE").is_none() {
            return;
        }
        let mut engine = GentleEngine::new();
        let result = engine
            .apply(Operation::FetchGenBankAccession {
                accession: "AY738222.1".to_string(),
                as_id: Some("pgl4_10_luc2".to_string()),
            })
            .expect("fetch and validate pGL4.10");
        let report = result
            .reporter_vector_validation
            .expect("fetch-time exact-vector validation");

        assert_eq!(report.status, ReporterVectorValidationStatus::Verified);
        assert_eq!(report.observed_length_bp, Some(4242));
        assert_eq!(report.observed_luc2_start_1based, Some(100));
        assert_eq!(
            report.observed_multiple_cloning_region.as_deref(),
            Some("1..70")
        );
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
        assert!(
            report.quarantined_records[0]
                .reasons
                .contains(&"unclear_or_unaccepted_license_status".to_string())
        );
    }

    #[test]
    fn reporter_annotation_reports_checksum_and_orf_sanity() {
        let mut record = test_record("r1", "ATGAAATAA");
        record.sequence_sha256 = Some("bad".to_string());
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
        assert!(
            report.quarantined_records[0]
                .reasons
                .contains(&"invalid_sequence_characters".to_string())
        );
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
        assert_eq!(result.biological_intent, "promoter_reporter_assay");
        assert_eq!(result.recommended_candidate_count, 1);
        assert_eq!(result.recommendations[0].reporter_id, "red");
        assert_eq!(result.rejected_candidate_count, 1);
        assert!(
            result.rejected_candidates[0]
                .reasons
                .iter()
                .any(|reason| reason.starts_with("color_mismatch"))
        );
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

    #[test]
    fn reporter_construct_handoff_uses_vkorc1_candidate_set_and_luciferase_route() {
        let engine = GentleEngine::new();
        let plan = engine
            .plan_reporter_construct_handoff(
                vkorc1_candidate_set_path(),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            )
            .expect("plan reporter construct handoff");
        assert_eq!(plan.schema, REPORTER_CONSTRUCT_HANDOFF_SCHEMA);
        assert_eq!(
            plan.biological_intent,
            "allele_paired_promoter_luciferase_reporter_handoff"
        );
        assert_eq!(
            plan.reporter_recommendation.biological_intent,
            "promoter_luciferase_reporter_assay"
        );
        assert_eq!(
            plan.provenance.macro_template_id,
            REPORTER_CONSTRUCT_MACRO_TEMPLATE_ID
        );
        assert_eq!(
            plan.selected_fragment.candidate_id,
            "vkorc1_rs9923231_context_enst00000498155_2413_3501_promoter_fragment"
        );
        assert_eq!(
            plan.selected_reporter
                .as_ref()
                .expect("selected reporter")
                .reporter_class,
            "luciferase"
        );
        assert!(plan.commands.iter().any(|command| {
            command
                .command
                .contains("allele_paired_promoter_luciferase_reporter")
        }));
    }

    #[test]
    fn reporter_construct_handoff_rejects_unknown_candidate_id() {
        let engine = GentleEngine::new();
        let err = engine
            .plan_reporter_construct_handoff(
                vkorc1_candidate_set_path(),
                Some("missing_candidate"),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            )
            .expect_err("unknown candidate should fail");
        assert_eq!(err.code, ErrorCode::InvalidInput);
    }

    #[test]
    fn reporter_construct_handoff_rejects_wrong_candidate_set_schema() {
        let td = tempdir().expect("tempdir");
        let path = td.path().join("candidate_set.json");
        std::fs::write(&path, r#"{"schema":"wrong.schema"}"#).expect("write candidate set");
        let engine = GentleEngine::new();
        let err = engine
            .plan_reporter_construct_handoff(
                &path.to_string_lossy(),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            )
            .expect_err("wrong schema should fail");
        assert_eq!(err.code, ErrorCode::InvalidInput);
    }

    #[test]
    fn reporter_construct_handoff_marks_non_luciferase_constraints_as_no_macro_route() {
        let constraints = ReporterConstraints {
            allowed_reporter_classes: vec!["fluorescent_protein".to_string()],
            ..ReporterConstraints::default()
        };
        let engine = GentleEngine::new();
        let plan = engine
            .plan_reporter_construct_handoff(
                vkorc1_candidate_set_path(),
                None,
                None,
                Some(constraints),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            )
            .expect("plan no compatible macro route");
        assert_eq!(plan.status, "no_compatible_macro_route");
        assert!(plan.selected_reporter.is_none());
        assert!(plan.port_bindings.is_empty());
        assert!(plan.commands.is_empty());
    }

    #[test]
    fn reporter_construct_handoff_reports_backbone_and_port_readiness() {
        let mut engine = GentleEngine::new();
        engine
            .apply(Operation::LoadFile {
                path: DEFAULT_REPORTER_BACKBONE_LOAD_PATH.to_string(),
                as_id: Some(DEFAULT_REPORTER_BACKBONE_SEQ_ID.to_string()),
            })
            .expect("load reporter backbone");
        let plan = engine
            .plan_reporter_construct_handoff(
                vkorc1_candidate_set_path(),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                Some("provided_reference"),
                None,
                None,
            )
            .expect("plan reporter construct handoff");
        assert_eq!(
            plan.backbone.status,
            ReporterBackboneResolutionStatus::ResolvedInState
        );
        let reference = plan
            .port_bindings
            .iter()
            .find(|binding| binding.port_id == "reference_fragment_seq_id")
            .expect("reference binding");
        assert_eq!(
            reference.status,
            PortBindingStatus::ProvidedMissingFromState
        );
        let alternate = plan
            .port_bindings
            .iter()
            .find(|binding| binding.port_id == "alternate_fragment_seq_id")
            .expect("alternate binding");
        assert_eq!(alternate.status, PortBindingStatus::Derivable);
        let backbone = plan
            .port_bindings
            .iter()
            .find(|binding| binding.port_id == "reporter_backbone_seq_id")
            .expect("backbone binding");
        assert_eq!(backbone.status, PortBindingStatus::Ready);
    }

    #[test]
    fn reporter_biological_intent_fields_are_additive_for_old_json() {
        let recommendation: ReporterRecommendationResult =
            serde_json::from_str(r#"{"schema":"gentle.reporter_recommendation.v1"}"#)
                .expect("decode old reporter recommendation");
        assert_eq!(recommendation.biological_intent, "");

        let handoff: ReporterConstructHandoffPlan =
            serde_json::from_str(r#"{"schema":"gentle.reporter_construct_handoff.v1"}"#)
                .expect("decode old reporter handoff");
        assert_eq!(handoff.biological_intent, "");
    }
}
