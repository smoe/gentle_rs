//! Provenance-preserving import and evaluation of externally supplied primer pairs.
//!
//! Source claims remain annotations. Transcript coverage, oligo QC, genomic
//! carryover and optional prepared-genome specificity are always recomputed by
//! GENtle's existing shared assay paths.

use std::{collections::BTreeMap, fs, path::Path};

use serde::Deserialize;

use super::*;

#[derive(Debug, Deserialize, Default)]
#[serde(default)]
struct ExternalPrimerPairTsvRow {
    source_kind: String,
    provider: String,
    #[serde(alias = "catalog_id")]
    catalogue_id: String,
    source_url: String,
    claimed_accession: String,
    aliases: String,
    forward_sequence_5_to_3: String,
    reverse_sequence_5_to_3: String,
    claimed_target: String,
    validation_claims: String,
    annotations_json: String,
}

#[derive(Debug)]
struct NormalizedExternalPrimerSource {
    forward: String,
    reverse: String,
    provenance: ExternalPrimerPairSourceProvenance,
}

impl ExternalPrimerPairSourceKind {
    fn as_external_import_str(self) -> &'static str {
        match self {
            Self::External => "external",
            Self::CommercialCatalogue => "commercial_catalogue",
            Self::Literature => "literature",
            Self::Laboratory => "laboratory",
        }
    }

    fn parse_external_import(raw: &str) -> Result<Self, EngineError> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "" | "external" | "imported_external" => Ok(Self::External),
            "commercial" | "commercial_catalog" | "commercial_catalogue" | "vendor" => {
                Ok(Self::CommercialCatalogue)
            }
            "literature" | "publication" => Ok(Self::Literature),
            "laboratory" | "lab" => Ok(Self::Laboratory),
            other => Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unknown external primer source_kind '{other}'; expected external, commercial_catalogue, literature, or laboratory"
                ),
                cause_chain: vec![],
            }),
        }
    }

    fn summary_origin(self) -> PrimerPairSummaryOrigin {
        match self {
            Self::CommercialCatalogue => PrimerPairSummaryOrigin::ImportedCommercial,
            Self::External | Self::Literature | Self::Laboratory => {
                PrimerPairSummaryOrigin::ImportedExternal
            }
        }
    }
}

impl GentleEngine {
    fn external_primer_split_list(raw: &str) -> Vec<String> {
        let mut values = raw
            .split(['|', ';'])
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(ToString::to_string)
            .collect::<Vec<_>>();
        values.sort();
        values.dedup();
        values
    }

    fn normalize_external_primer_sequence(
        raw: &str,
        row_number: usize,
        role: &str,
    ) -> Result<String, EngineError> {
        let mut normalized = String::new();
        for (index, character) in raw.chars().enumerate() {
            if character.is_whitespace() || character.is_ascii_digit() {
                continue;
            }
            if !character.is_ascii() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "External primer row {row_number} {role} sequence contains non-ASCII character '{}' at input character {}",
                        character,
                        index.saturating_add(1)
                    ),
                    cause_chain: vec![],
                });
            }
            let byte = character as u8;
            if !IupacCode::is_valid_letter(byte) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "External primer row {row_number} {role} sequence contains invalid IUPAC character '{}' at input character {}",
                        character,
                        index.saturating_add(1)
                    ),
                    cause_chain: vec![],
                });
            }
            let upper = byte.to_ascii_uppercase();
            normalized.push(if upper == b'U' { 'T' } else { upper as char });
        }
        if normalized.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "External primer row {row_number} {role} sequence contains no IUPAC bases after normalization"
                ),
                cause_chain: vec![],
            });
        }
        Ok(normalized)
    }

    fn external_primer_batch_from_tsv(
        bytes: &[u8],
        default_batch_id: &str,
    ) -> Result<ExternalPrimerPairBatch, EngineError> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(true)
            .from_reader(bytes);
        let mut pairs = vec![];
        for (index, row) in reader.deserialize::<ExternalPrimerPairTsvRow>().enumerate() {
            let row = row.map_err(|error| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not parse external primer TSV data row {}: {error}",
                    index.saturating_add(1)
                ),
                cause_chain: vec![],
            })?;
            let annotations = if row.annotations_json.trim().is_empty() {
                BTreeMap::new()
            } else {
                serde_json::from_str::<BTreeMap<String, String>>(&row.annotations_json).map_err(
                    |error| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "External primer TSV data row {} has invalid annotations_json: {error}",
                            index.saturating_add(1)
                        ),
                        cause_chain: vec![],
                    },
                )?
            };
            pairs.push(ExternalPrimerPairInput {
                source_kind: ExternalPrimerPairSourceKind::parse_external_import(&row.source_kind)?,
                provider: row.provider,
                catalogue_id: row.catalogue_id,
                source_url: row.source_url,
                claimed_accession: row.claimed_accession,
                aliases: Self::external_primer_split_list(&row.aliases),
                forward_sequence_5_to_3: row.forward_sequence_5_to_3,
                reverse_sequence_5_to_3: row.reverse_sequence_5_to_3,
                claimed_target: row.claimed_target,
                validation_claims: Self::external_primer_split_list(&row.validation_claims),
                annotations,
            });
        }
        Ok(ExternalPrimerPairBatch {
            schema: EXTERNAL_PRIMER_PAIR_BATCH_SCHEMA.to_string(),
            batch_id: default_batch_id.to_string(),
            pairs,
        })
    }

    /// Load a JSON or TSV external-primer batch and attach file provenance.
    pub fn load_external_primer_pair_batch(
        path: &str,
        requested_format: Option<&str>,
    ) -> Result<(ExternalPrimerPairBatch, ExternalPrimerPairBatchProvenance), EngineError> {
        let bytes = fs::read(path).map_err(|error| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read external primer-pair batch '{path}': {error}"),
            cause_chain: vec![],
        })?;
        let inferred = Path::new(path)
            .extension()
            .and_then(|value| value.to_str())
            .unwrap_or_default()
            .to_ascii_lowercase();
        let format = requested_format
            .map(str::trim)
            .filter(|value| !value.is_empty() && !value.eq_ignore_ascii_case("auto"))
            .map(|value| value.to_ascii_lowercase())
            .unwrap_or_else(|| {
                if matches!(inferred.as_str(), "tsv" | "tab") {
                    "tsv".to_string()
                } else {
                    "json".to_string()
                }
            });
        let source_sha256 = sha256_prefixed_bytes(&bytes);
        let default_batch_id = short_sha256_id("external_primer_batch", &source_sha256);
        let mut batch = match format.as_str() {
            "json" => {
                serde_json::from_slice::<ExternalPrimerPairBatch>(&bytes).map_err(|error| {
                    EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not parse external primer-pair JSON batch '{path}': {error}"
                        ),
                        cause_chain: vec![],
                    }
                })?
            }
            "tsv" => Self::external_primer_batch_from_tsv(&bytes, &default_batch_id)?,
            other => {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Unknown external primer-pair format '{other}'; expected auto, json, or tsv"
                    ),
                    cause_chain: vec![],
                });
            }
        };
        if batch.schema.trim().is_empty() {
            batch.schema = EXTERNAL_PRIMER_PAIR_BATCH_SCHEMA.to_string();
        } else if batch.schema != EXTERNAL_PRIMER_PAIR_BATCH_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "External primer-pair batch schema '{}' is unsupported; expected '{}'",
                    batch.schema, EXTERNAL_PRIMER_PAIR_BATCH_SCHEMA
                ),
                cause_chain: vec![],
            });
        }
        if batch.batch_id.trim().is_empty() {
            batch.batch_id = default_batch_id;
        }
        if batch.pairs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "External primer-pair batch contains no pair rows".to_string(),
                cause_chain: vec![],
            });
        }
        Ok((
            batch,
            ExternalPrimerPairBatchProvenance {
                input_format: format,
                source_path: path.to_string(),
                source_sha256,
            },
        ))
    }

    fn normalize_external_primer_sources(
        request: &ExternalPrimerPairImportRequest,
    ) -> Result<Vec<NormalizedExternalPrimerSource>, EngineError> {
        if request.batch.schema != EXTERNAL_PRIMER_PAIR_BATCH_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "External primer-pair batch schema '{}' is unsupported; expected '{}'",
                    request.batch.schema, EXTERNAL_PRIMER_PAIR_BATCH_SCHEMA
                ),
                cause_chain: vec![],
            });
        }
        if request.batch.pairs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "External primer-pair batch contains no pair rows".to_string(),
                cause_chain: vec![],
            });
        }
        let mut normalized = Vec::with_capacity(request.batch.pairs.len());
        for (index, input) in request.batch.pairs.iter().enumerate() {
            let row_number = index.saturating_add(1);
            if input.provider.trim().is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "External primer row {row_number} requires a non-empty provider"
                    ),
                    cause_chain: vec![],
                });
            }
            if input.source_kind == ExternalPrimerPairSourceKind::CommercialCatalogue
                && input.catalogue_id.trim().is_empty()
            {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "External primer row {row_number} declares commercial_catalogue origin but has no catalogue_id"
                    ),
                    cause_chain: vec![],
                });
            }
            let forward = Self::normalize_external_primer_sequence(
                &input.forward_sequence_5_to_3,
                row_number,
                "forward",
            )?;
            let reverse = Self::normalize_external_primer_sequence(
                &input.reverse_sequence_5_to_3,
                row_number,
                "reverse",
            )?;
            let mut aliases = input
                .aliases
                .iter()
                .map(|value| value.trim())
                .filter(|value| !value.is_empty())
                .map(ToString::to_string)
                .collect::<Vec<_>>();
            aliases.sort();
            aliases.dedup();
            let mut validation_claims = input
                .validation_claims
                .iter()
                .map(|value| value.trim())
                .filter(|value| !value.is_empty())
                .map(ToString::to_string)
                .collect::<Vec<_>>();
            validation_claims.sort();
            validation_claims.dedup();
            let source_key = serde_json::to_string(&(
                input.source_kind.as_external_import_str(),
                input.provider.trim(),
                input.catalogue_id.trim(),
                input.source_url.trim(),
                input.claimed_accession.trim(),
                &aliases,
                input.claimed_target.trim(),
                &validation_claims,
                &input.annotations,
                &forward,
                &reverse,
            ))
            .map_err(|error| EngineError {
                code: ErrorCode::Internal,
                message: format!(
                    "Could not derive a stable identity for external primer row {row_number}: {error}"
                ),
                cause_chain: vec![],
            })?;
            normalized.push(NormalizedExternalPrimerSource {
                forward,
                reverse,
                provenance: ExternalPrimerPairSourceProvenance {
                    source_record_id: short_sha256_id("external_primer_source", &source_key),
                    input_row_number: row_number,
                    source_kind: input.source_kind,
                    provider: input.provider.trim().to_string(),
                    catalogue_id: input.catalogue_id.trim().to_string(),
                    source_url: input.source_url.trim().to_string(),
                    claimed_accession: input.claimed_accession.trim().to_string(),
                    aliases,
                    claimed_target: input.claimed_target.trim().to_string(),
                    validation_claims,
                    annotations: input.annotations.clone(),
                    input_format: request.input_provenance.input_format.clone(),
                    source_path: request.input_provenance.source_path.clone(),
                    source_sha256: request.input_provenance.source_sha256.clone(),
                    claim_evidence_status: "provenance_only_not_used_for_coverage_or_specificity"
                        .to_string(),
                },
            });
        }
        Ok(normalized)
    }

    fn normalized_external_primer_batch_sha256(
        sources: &[NormalizedExternalPrimerSource],
    ) -> Result<String, EngineError> {
        let mut rows = sources
            .iter()
            .map(|source| {
                (
                    source.provenance.source_record_id.as_str(),
                    source.forward.as_str(),
                    source.reverse.as_str(),
                )
            })
            .collect::<Vec<_>>();
        rows.sort_unstable();
        let bytes = serde_json::to_vec(&rows).map_err(|error| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize normalized external primer batch identity: {error}"
            ),
            cause_chain: vec![],
        })?;
        Ok(sha256_prefixed_bytes(&bytes))
    }

    fn external_primer_import_report_id(
        request: &ExternalPrimerPairImportRequest,
        normalized_batch_sha256: &str,
    ) -> Result<String, EngineError> {
        let identity = serde_json::to_vec(&(
            normalized_batch_sha256,
            request.seq_id.as_str(),
            request.source_feature_id,
            request.transcript_id.as_deref(),
            request.min_amplicon_bp,
            request.max_amplicon_bp,
            request.max_mismatches,
            request.require_3prime_exact_bases,
            request.transcript_order,
            request.transcript_map_coordinate_mode,
            request.specificity.as_ref(),
            request.materialize_products,
            &request.product_gel_ladders,
        ))
        .map_err(|error| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize external primer import evaluation identity: {error}"
            ),
            cause_chain: vec![],
        })?;
        Ok(short_sha256_id(
            "external_primer_import",
            &sha256_prefixed_bytes(&identity),
        ))
    }

    fn external_primer_oligo_assessment(
        role: &str,
        sequence: &str,
        origin: PrimerPairSummaryOrigin,
        cdna_anneal_hit_count: usize,
        qc: &OligoQcReport,
    ) -> ExternalPrimerOligoAssessment {
        let qc_role = format!("{role}_primer");
        let qc_row = qc.oligos.iter().find(|row| row.role == qc_role);
        let metrics = Self::compute_primer_heuristic_metrics(sequence.as_bytes());
        let gc_fraction = Self::sequence_gc_fraction(sequence.as_bytes()).unwrap_or(0.0);
        ExternalPrimerOligoAssessment {
            primer_id: short_sha256_id("primer", sequence),
            role: role.to_string(),
            origin,
            sequence_5_to_3: sequence.to_string(),
            length_nt: sequence.len(),
            anneal_length_nt: sequence.len(),
            tm_c: Self::estimate_primer_tm_c(sequence.as_bytes()),
            tm_method: "gentle.shared_primer_tm.v1".to_string(),
            tm_assumptions: Self::primer_tm_model_description(),
            gc_fraction,
            gc_percent: gc_fraction * 100.0,
            cdna_anneal_hit_count,
            three_prime_base: char::from(metrics.three_prime_base).to_string(),
            three_prime_gc_clamp: metrics.three_prime_gc_clamp,
            longest_homopolymer_run_bp: metrics.longest_homopolymer_run_bp,
            self_complementary_run_bp: metrics.self_complementary_run_bp,
            self_3prime_complementary_run_bp: qc_row
                .map(|row| row.self_3prime_complementary_run_bp)
                .unwrap_or_default(),
            qc_status: qc_row
                .map(|row| row.status.clone())
                .unwrap_or_else(|| "unknown".to_string()),
            qc_warnings: qc_row.map(|row| row.warnings.clone()).unwrap_or_default(),
        }
    }

    fn external_primer_specificity_assessment(
        &self,
        forward: &str,
        reverse: &str,
        request: Option<&ExternalPrimerPairSpecificityRequest>,
    ) -> ExternalPrimerPairSpecificityAssessment {
        let Some(request) = request else {
            return ExternalPrimerPairSpecificityAssessment {
                status: "not_run".to_string(),
                reason: "No prepared-genome specificity request was supplied; vendor claims are not treated as a pass."
                    .to_string(),
                report: None,
            };
        };
        match self.assess_primer_pair_specificity(
            None,
            None,
            None,
            Some(forward),
            Some(reverse),
            &request.target_genome_id,
            request.policy.clone(),
            request.catalog_path.as_deref(),
            request.cache_dir.as_deref(),
        ) {
            Ok(report) => ExternalPrimerPairSpecificityAssessment {
                status: if !report.search_completeness.complete {
                    "incomplete"
                } else if report.summary.specificity_pass {
                    "pass"
                } else {
                    "specificity_fail"
                }
                .to_string(),
                reason: report.summary.summary.clone(),
                report: Some(Box::new(report)),
            },
            Err(error) => ExternalPrimerPairSpecificityAssessment {
                status: "error".to_string(),
                reason: error.to_string(),
                report: None,
            },
        }
    }

    fn external_primer_artifact_path(directory: &str, pair_id: &str, suffix: &str) -> String {
        Path::new(directory)
            .join(format!("{pair_id}.{suffix}"))
            .to_string_lossy()
            .to_string()
    }

    fn write_external_primer_json<T: Serialize>(
        path: &str,
        value: &T,
        label: &str,
    ) -> Result<(), EngineError> {
        if let Some(parent) = Path::new(path).parent()
            && !parent.as_os_str().is_empty()
        {
            fs::create_dir_all(parent).map_err(|error| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not create {label} directory '{}': {error}",
                    parent.display()
                ),
                cause_chain: vec![],
            })?;
        }
        let bytes = serde_json::to_vec_pretty(value).map_err(|error| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize {label}: {error}"),
            cause_chain: vec![],
        })?;
        fs::write(path, bytes).map_err(|error| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write {label} '{path}': {error}"),
            cause_chain: vec![],
        })
    }

    /// Import, deduplicate and evaluate one batch of externally supplied pairs.
    pub(super) fn import_external_primer_pairs(
        &mut self,
        op_result: &mut OpResult,
        mut request: ExternalPrimerPairImportRequest,
    ) -> Result<ExternalPrimerPairImportReport, EngineError> {
        let normalized_sources = Self::normalize_external_primer_sources(&request)?;
        let normalized_batch_sha256 =
            Self::normalized_external_primer_batch_sha256(&normalized_sources)?;
        if let Some(specificity) = request.specificity.as_ref()
            && specificity.target_genome_id.trim().is_empty()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "External primer specificity target_genome_id cannot be empty".to_string(),
                cause_chain: vec![],
            });
        }
        if request.report_id.trim().is_empty() {
            request.report_id =
                Self::external_primer_import_report_id(&request, &normalized_batch_sha256)?;
        }
        request.report_id = Self::normalize_primer_design_report_id(&request.report_id)?;

        let mut grouped = BTreeMap::<String, Vec<NormalizedExternalPrimerSource>>::new();
        for source in normalized_sources {
            let pair_id = short_sha256_id(
                "primer_pair",
                &format!("{}\n{}", source.forward, source.reverse),
            );
            grouped.entry(pair_id).or_default().push(source);
        }
        let source_record_count = request.batch.pairs.len();
        let duplicate_source_record_count = source_record_count.saturating_sub(grouped.len());
        let mut pair_results = Vec::with_capacity(grouped.len());
        let artifact_directory = request
            .artifact_output_dir
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(ToString::to_string);
        if let Some(directory) = artifact_directory.as_deref() {
            fs::create_dir_all(directory).map_err(|error| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not create external primer artifact directory '{directory}': {error}"
                ),
                cause_chain: vec![],
            })?;
        }

        for (pair_id, mut sources) in grouped {
            sources.sort_by(|left, right| {
                left.provenance
                    .source_record_id
                    .cmp(&right.provenance.source_record_id)
            });
            let forward = sources[0].forward.clone();
            let reverse = sources[0].reverse.clone();
            let source_rows = sources
                .iter()
                .map(|source| source.provenance.clone())
                .collect::<Vec<_>>();
            let mut aliases = source_rows
                .iter()
                .flat_map(|source| source.aliases.iter().cloned())
                .collect::<Vec<_>>();
            aliases.sort();
            aliases.dedup();
            let source_origins = source_rows
                .iter()
                .map(|source| source.source_kind.summary_origin())
                .collect::<Vec<_>>();
            let mut origins = vec![];
            if source_origins
                .iter()
                .any(|origin| *origin == PrimerPairSummaryOrigin::ImportedCommercial)
            {
                origins.push(PrimerPairSummaryOrigin::ImportedCommercial);
            }
            if source_origins
                .iter()
                .any(|origin| *origin == PrimerPairSummaryOrigin::ImportedExternal)
            {
                origins.push(PrimerPairSummaryOrigin::ImportedExternal);
            }
            let oligo_origin = if origins.contains(&PrimerPairSummaryOrigin::ImportedCommercial) {
                PrimerPairSummaryOrigin::ImportedCommercial
            } else {
                PrimerPairSummaryOrigin::ImportedExternal
            };

            let cdna_assay = self.test_cdna_pcr_assay_with_map_options(
                &request.seq_id,
                request.source_feature_id,
                &forward,
                &reverse,
                request.transcript_id.as_deref(),
                request.min_amplicon_bp,
                request.max_amplicon_bp,
                request.max_mismatches,
                request.require_3prime_exact_bases,
                request.transcript_order,
                request.transcript_map_coordinate_mode,
            )?;
            let forward_hit_count = cdna_assay
                .transcript_results
                .iter()
                .map(|row| row.forward_hits.len())
                .sum();
            let reverse_hit_count = cdna_assay
                .transcript_results
                .iter()
                .map(|row| row.reverse_hits.len())
                .sum();
            let oligo_qc = Self::build_oligo_qc_report(
                "pcr",
                &[
                    ("forward_primer", "forward_primer", forward.as_str()),
                    ("reverse_primer", "reverse_primer", reverse.as_str()),
                ],
            );
            let forward_assessment = Self::external_primer_oligo_assessment(
                "forward",
                &forward,
                oligo_origin,
                forward_hit_count,
                &oligo_qc,
            );
            let reverse_assessment = Self::external_primer_oligo_assessment(
                "reverse",
                &reverse,
                oligo_origin,
                reverse_hit_count,
                &oligo_qc,
            );
            let specificity = self.external_primer_specificity_assessment(
                &forward,
                &reverse,
                request.specificity.as_ref(),
            );
            let mut artifacts = ExternalPrimerPairArtifacts::default();
            if let Some(directory) = artifact_directory.as_deref() {
                let report_path =
                    Self::external_primer_artifact_path(directory, &pair_id, "cdna_assay.json");
                Self::write_external_primer_json(
                    &report_path,
                    &cdna_assay,
                    "external primer cDNA assay report",
                )?;
                artifacts.cdna_report_json_path = Some(report_path);
                let map_path =
                    Self::external_primer_artifact_path(directory, &pair_id, "transcript_map.svg");
                let map_svg = cdna_assay
                    .transcript_map
                    .as_ref()
                    .map(|map| map.svg.as_str())
                    .unwrap_or_default();
                fs::write(&map_path, map_svg).map_err(|error| EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "Could not write external primer transcript map '{map_path}': {error}"
                    ),
                    cause_chain: vec![],
                })?;
                artifacts.transcript_map_svg_path = Some(map_path);
            }
            let product_materialization = if request.materialize_products {
                let gel_path = artifact_directory.as_deref().map(|directory| {
                    Self::external_primer_artifact_path(directory, &pair_id, "product_gel.svg")
                });
                let prefix = format!("external_{pair_id}");
                let materialization = self.materialize_cdna_assay_products(
                    op_result,
                    &cdna_assay,
                    Some(&prefix),
                    gel_path.as_deref(),
                    (!request.product_gel_ladders.is_empty())
                        .then_some(request.product_gel_ladders.as_slice()),
                )?;
                artifacts.product_gel_svg_path = materialization.product_gel_svg_path.clone();
                Some(materialization)
            } else {
                None
            };
            let mut warnings = cdna_assay.warnings.clone();
            if specificity.status == "error" {
                warnings.push(format!(
                    "Whole-genome specificity failed to execute: {}",
                    specificity.reason
                ));
            }
            pair_results.push(ExternalPrimerPairAssessment {
                pair_id,
                aliases,
                origins,
                duplicate_source_record_count: source_rows.len().saturating_sub(1),
                sources: source_rows,
                tm_delta_c: (forward_assessment.tm_c - reverse_assessment.tm_c).abs(),
                forward: forward_assessment,
                reverse: reverse_assessment,
                oligo_qc,
                cdna_assay,
                specificity,
                artifacts,
                product_materialization,
                vendor_claims_used_as_biological_evidence: false,
                warnings,
            });
        }

        let mut warnings = vec![];
        if request.specificity.is_none() {
            warnings.push(
                "Whole-genome specificity was not run; imported source claims were retained as provenance and do not imply a pass."
                    .to_string(),
            );
        }
        if request.materialize_products && artifact_directory.is_none() {
            warnings.push(
                "Products were materialized, but no artifact_output_dir was supplied, so product gels were not written."
                    .to_string(),
            );
        }
        let report = ExternalPrimerPairImportReport {
            schema: EXTERNAL_PRIMER_PAIR_IMPORT_REPORT_SCHEMA.to_string(),
            report_id: request.report_id.clone(),
            batch_id: request.batch.batch_id.clone(),
            normalized_batch_sha256,
            generated_at_unix_ms: Self::now_unix_ms(),
            seq_id: request.seq_id,
            source_feature_id: request.source_feature_id,
            input_provenance: request.input_provenance,
            source_record_count,
            unique_pair_count: pair_results.len(),
            duplicate_source_record_count,
            pairs: pair_results,
            warnings,
        };
        let mut store = self.read_primer_design_store();
        store
            .external_primer_pair_imports
            .insert(report.report_id.clone(), report.clone());
        self.write_primer_design_store(store)?;
        Ok(report)
    }

    pub fn get_external_primer_pair_import_report(
        &self,
        report_id: &str,
    ) -> Result<ExternalPrimerPairImportReport, EngineError> {
        let report_id = Self::normalize_primer_design_report_id(report_id)?;
        self.read_primer_design_store()
            .external_primer_pair_imports
            .get(&report_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("External primer-pair import report '{report_id}' not found"),
                cause_chain: vec![],
            })
    }

    pub fn list_external_primer_pair_import_report_ids(&self) -> Vec<String> {
        let mut report_ids = self
            .read_primer_design_store()
            .external_primer_pair_imports
            .keys()
            .cloned()
            .collect::<Vec<_>>();
        report_ids.sort();
        report_ids
    }

    pub fn export_external_primer_pair_import_report(
        &self,
        report_id: &str,
        path: &str,
    ) -> Result<ExternalPrimerPairImportReport, EngineError> {
        let report = self.get_external_primer_pair_import_report(report_id)?;
        Self::write_external_primer_json(path, &report, "external primer-pair import report")?;
        Ok(report)
    }
}
