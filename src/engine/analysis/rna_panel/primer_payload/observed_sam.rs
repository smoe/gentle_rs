//! S2.3 verifies bound SAM/annotation evidence and publishes separate, immutable reports.

use super::*;
use crate::observed_transcript_structures as kernel;
use serde_json::json;

#[path = "observed_sam/publication.rs"]
mod publication;
#[cfg(test)]
#[path = "observed_sam/tests.rs"]
mod tests;

const CAP: u64 = 128 * 1024 * 1024;
const TABLES: [&str; 6] = [
    "query_structures.tsv",
    "alignment_observations.tsv",
    "structure_catalog.tsv",
    "reference_compatibility.tsv",
    "junction_catalog.tsv",
    "unavailable.tsv",
];
const DATA: &str = "analysis.json";

type Request = SummarizeObservedTranscriptStructuresFromSamRequest;

#[cfg(test)]
impl GentleEngine {
    pub(crate) fn observed_sam_test_fixture() -> (tempfile::TempDir, Request) {
        tests::fixture()
    }
}

fn digest<T: Serialize>(value: &T) -> Result<String, EngineError> {
    crate::digest_utils::canonical_json_string(value)
        .map(|s| sha256_prefixed_bytes(s.as_bytes()))
        .map_err(|e| bad(e.to_string()))
}

struct Size {
    bytes: u64,
    limit: u64,
}
impl Write for Size {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len() as u64)
            .filter(|n| *n <= self.limit)
            .ok_or_else(|| std::io::Error::other("observed_sam_json_memory_budget"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn size<T: Serialize>(value: &T, limit: u64) -> Result<u64, EngineError> {
    let mut counter = Size { bytes: 0, limit };
    serde_json::to_writer(&mut counter, value).map_err(|e| bad(e.to_string()))?;
    Ok(counter.bytes)
}

fn policy(r: &Request) -> Result<String, EngineError> {
    digest(&(
        "gentle.observed_sam_policy.v1",
        "masked_supported_bases_all_best_equivalent_paths",
        "exact_physical_junctions",
        r.version_policy,
        r.unpadded_reference,
        r.required_flank_bases,
        &r.text_budgets,
        &r.kernel_budgets,
        &r.budgets,
    ))
}

fn input_files(r: &Request) -> Vec<&PrimerPayloadFile> {
    let mut files = vec![&r.sam, &r.manifest, &r.dictionary];
    for source in &r.annotation_sources {
        match source {
            ObservedSamAnnotationBinding::Prepared {
                report, annotation, ..
            } => {
                files.push(report);
                files.extend(annotation.iter());
            }
            ObservedSamAnnotationBinding::Normalized { annotation, .. } => files.push(annotation),
        }
    }
    files
}

fn absolute(raw: &str) -> Result<PathBuf, EngineError> {
    let p = kernel::validate_sam_raw_path(raw)?;
    Ok(if p.is_absolute() {
        p.to_owned()
    } else {
        std::env::current_dir()
            .map_err(|e| bad(e.to_string()))?
            .join(p)
    })
}

fn paths(r: &Request) -> Result<PathBuf, EngineError> {
    let out = absolute(&r.output_dir)?;
    if let Some(marker) = &r.cancel_path {
        let marker = absolute(marker)?;
        let folded = |p: &Path| p.to_string_lossy().to_lowercase();
        if marker == out
            || marker.starts_with(&out)
            || ((cfg!(windows) || cfg!(target_os = "macos"))
                && (folded(&marker) == folded(&out)
                    || folded(&marker)
                        .starts_with(&(folded(&out) + std::path::MAIN_SEPARATOR_STR))))
        {
            return Err(bad("observed_sam_control_collision"));
        }
        for file in input_files(r) {
            let input = absolute(&file.path)?;
            if input == marker
                || same_file::is_same_file(&input, &marker).unwrap_or(false)
                || ((cfg!(windows) || cfg!(target_os = "macos"))
                    && folded(&input) == folded(&marker))
            {
                return Err(bad("observed_sam_control_collision"));
            }
        }
    }
    Ok(out)
}

fn admit(r: &Request) -> Result<(), EngineError> {
    let b = &r.budgets;
    r.kernel_budgets.validate(r.required_flank_bases)?;
    if r.schema != OBSERVED_SAM_REQUEST_SCHEMA
        || !(1..=CAP).contains(&b.max_input_bytes)
        || !(1..=CAP).contains(&b.max_annotation_bytes)
        || !(1..=CAP).contains(&b.max_output_bytes)
        || !(1..=8 * 1024 * 1024).contains(&b.max_output_row_bytes)
        || !(1..=CAP as usize).contains(&b.max_retained_bytes)
        || b.max_retained_bytes <= 65536
        || !(1..=4096).contains(&b.max_units)
        || !(1..=100_000_000).contains(&b.max_projection_work)
        || r.text_budgets.max_retained_bytes > b.max_retained_bytes
        || r.text_budgets.max_input_bytes as u64 > b.max_input_bytes
        || r.annotation_sources.len() > 1000
    {
        return Err(bad("observed_sam_admission"));
    }
    let mut total = 0u64;
    for file in input_files(r) {
        kernel::validate_sam_raw_path(&file.path)?;
        if !is_hash(&file.sha256) {
            return Err(bad("observed_sam_input_binding"));
        }
        total = total
            .checked_add(file.bytes)
            .filter(|n| *n <= b.max_input_bytes)
            .ok_or_else(|| bad("observed_sam_input_budget"))?;
    }
    kernel::validate_sam_raw_path(&r.output_dir)?;
    if let Some(marker) = &r.cancel_path {
        kernel::validate_sam_raw_path(marker)?;
    }
    let selectors = std::iter::once(&r.reference_gene)
        .chain(&r.comparison_genes)
        .collect::<Vec<_>>();
    let mut expected = BTreeMap::new();
    for selector in selectors {
        if expected
            .insert(&selector.source_id, &selector.gene_id)
            .is_some()
        {
            return Err(bad("observed_sam_selector_binding_mismatch"));
        }
    }
    let mut actual = BTreeMap::new();
    for source in &r.annotation_sources {
        let selector = source.selector();
        if actual
            .insert(&selector.source_id, &selector.gene_id)
            .is_some()
            || expected.get(&selector.source_id) != Some(&&selector.gene_id)
        {
            return Err(bad("observed_sam_selector_binding_mismatch"));
        }
    }
    if expected != actual {
        return Err(bad("observed_sam_selector_binding_mismatch"));
    }
    size(r, (b.max_retained_bytes as u64 - 65536) / 64)?;
    Ok(())
}

fn json_bound<T: serde::de::DeserializeOwned>(
    file: &PrimerPayloadFile,
    r: &Request,
    ctl: &mut Control<'_>,
    retained: &mut usize,
) -> Result<T, EngineError> {
    if file.bytes > r.budgets.max_annotation_bytes {
        return Err(bad("observed_sam_annotation_budget"));
    }
    let charge = usize::try_from(file.bytes)
        .ok()
        .and_then(|n| n.checked_mul(64))
        .and_then(|n| n.checked_add(65536))
        .ok_or_else(|| bad("observed_sam_annotation_memory_budget"))?;
    *retained = retained
        .checked_add(charge)
        .filter(|n| *n < r.budgets.max_retained_bytes)
        .ok_or_else(|| bad("observed_sam_annotation_memory_budget"))?;
    kernel::read_sam_bound_json(file, r.budgets.max_retained_bytes, &mut || ctl.check())
}

fn annotations(
    r: &Request,
    ctl: &mut Control<'_>,
) -> Result<
    (
        Vec<kernel::AnnotationSource>,
        Vec<ObservedSamAnnotationEvidence>,
        BTreeMap<String, ObservedSamAnnotationAvailability>,
        usize,
    ),
    EngineError,
> {
    let mut sources = Vec::new();
    let mut evidence = Vec::new();
    let mut states = BTreeMap::new();
    let mut retained = 0usize;
    for binding in &r.annotation_sources {
        ctl.check()?;
        let selector = binding.selector();
        let (availability, preparation, annotation) = match binding {
            ObservedSamAnnotationBinding::Normalized {
                annotation,
                provider,
                ..
            } => {
                let value: PrimerPayloadAnnotation = json_bound(annotation, r, ctl, &mut retained)?;
                (
                    ObservedSamAnnotationAvailability::Available,
                    None,
                    Some((value, annotation.sha256.clone(), provider.clone())),
                )
            }
            ObservedSamAnnotationBinding::Prepared {
                report, annotation, ..
            } => {
                let prep: ObservedTranscriptAnnotationReport =
                    json_bound(report, r, ctl, &mut retained)?;
                annotation_prepare::validate_receipt_metadata(&prep)?;
                if encode(&prep.request.selector)? != encode(selector)?
                    || encode(&prep.request.reference)? != encode(&r.reference)?
                    || prep.request.reference_id != r.reference_id
                {
                    return Err(bad("observed_sam_selector_binding_mismatch"));
                }
                let usable = prep.status == ObservedAnnotationExecutionStatus::Completed
                    && prep
                        .analysis
                        .as_ref()
                        .is_some_and(|a| a.outcome == ObservedAnnotationOutcome::Completed);
                if usable {
                    let file = annotation
                        .as_ref()
                        .ok_or_else(|| bad("observed_sam_annotation_binding"))?;
                    let artifact = prep
                        .annotation
                        .as_ref()
                        .ok_or_else(|| bad("observed_sam_annotation_binding"))?;
                    if file.sha256 != artifact.sha256 || file.bytes != artifact.bytes {
                        return Err(bad("observed_sam_annotation_binding"));
                    }
                    let value: PrimerPayloadAnnotation = json_bound(file, r, ctl, &mut retained)?;
                    if value.release != format!("snapshot-{}", prep.request.gtf.sha256) {
                        return Err(bad("observed_sam_annotation_binding"));
                    }
                    let provider = prep.request.source.provider.clone();
                    (
                        ObservedSamAnnotationAvailability::Available,
                        Some(prep),
                        Some((value, file.sha256.clone(), provider)),
                    )
                } else {
                    if annotation.is_some() {
                        return Err(bad("observed_sam_unavailable_annotation_binding"));
                    }
                    let state = if prep.status == ObservedAnnotationExecutionStatus::Completed
                        && prep
                            .analysis
                            .as_ref()
                            .is_some_and(|a| a.outcome == ObservedAnnotationOutcome::MissingGene)
                    {
                        ObservedSamAnnotationAvailability::MissingGene
                    } else {
                        ObservedSamAnnotationAvailability::Refused
                    };
                    (state, Some(prep), None)
                }
            }
        };
        if let Some((annotation, sha256, provider)) = annotation {
            if annotation.schema != PRIMER_PAYLOAD_ANNOTATION_SCHEMA
                || annotation.reference_id != r.reference_id
                || annotation.assembly != r.reference.assembly
                || annotation.reference_fasta_sha256 != r.reference.fasta_sha256
                || annotation.coordinate_convention != "zero_based_half_open"
                || annotation.transcripts.is_empty()
            {
                return Err(bad("observed_sam_annotation_binding"));
            }
            sources.push(kernel::AnnotationSource {
                source_id: selector.source_id.clone(),
                provider,
                annotation_sha256: sha256,
                annotation,
            });
        }
        states.insert(selector.source_id.clone(), availability);
        evidence.push(ObservedSamAnnotationEvidence {
            binding: binding.clone(),
            availability,
            preparation,
        });
    }
    Ok((sources, evidence, states, retained))
}

fn report_counts(data: &ObservedSamAnalysis) -> ObservedSamCounts {
    let assessed = data
        .units
        .iter()
        .filter(|u| u.all_query_views_assessed)
        .count();
    ObservedSamCounts {
        query_units: data.units.len(),
        assessed_units: assessed,
        unavailable_units: data.units.len() - assessed,
        observations: data.observations.len(),
        optimal_observations: data
            .observations
            .iter()
            .filter(|o| o.on_optimal_path)
            .count(),
        chain_families: data.chain_families.len(),
        junctions: data.junctions.len(),
    }
}

fn nonclaims() -> Vec<String> {
    ["Counts describe producer-declared SAM query units, not verified original reads or independent molecules.",
     "EOF, retained artifacts and hashes do not verify exhaustive aligner search, the original BAM or FASTA bytes.",
     "Masks, selection scope, conversion and search coverage retain their explicit producer provenance.",
     "Mapping orientation is not RNA strand; junction geometry is not novelty, expression, mechanism or protein identity.",
     "Only retained best-equivalent paths contribute to catalogs; longer unspliced competitors may displace spliced mappings.",
     "Unknown primer masks do not become primer-free qualified evidence; qualification and mapping certainty are independent."]
        .into_iter().map(str::to_owned).collect()
}

impl GentleEngine {
    /// Verify bound SAM/query units and explicit annotations; no sequence or project mutation.
    pub fn summarize_observed_transcript_structures_from_sam(
        &self,
        r: Request,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<ObservedSamReport, EngineError> {
        if matches!(
            &r.conversion,
            kernel::ObservedSamConversion::ExecutedIndexedSelection { .. }
        ) {
            return Err(bad(
                "observed_sam_executed_selection_requires_local_execution",
            ));
        }
        self.summarize_observed_sam_internal(r, progress)
    }

    pub(crate) fn summarize_observed_sam_internal(
        &self,
        r: Request,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<ObservedSamReport, EngineError> {
        admit(&r)?;
        let out = paths(&r)?;
        if out.parent().is_none_or(|p| !p.is_dir()) {
            return Err(bad("observed_sam_output_parent"));
        }
        if fs::symlink_metadata(&out).is_ok() {
            return Err(bad("observed_sam_output_exists"));
        }
        let identity = digest(&r)?;
        let mut ctl = Control::new(r.cancel_path.clone(), identity.clone(), 5, progress)?;
        let mut data = ObservedSamData {
            schema: OBSERVED_SAM_DATA_SCHEMA.into(),
            request_sha256: identity.clone(),
            intake: None,
            annotations: Vec::new(),
            analysis: None,
        };
        let result = (|| {
            let (sources, evidence, states, retained) = annotations(&r, &mut ctl)?;
            data.annotations = evidence;
            ctl.completed = 1;
            if ctl.stop(true) {
                return Err(bad(INTERRUPTED));
            }
            let mut limits = r.text_budgets.clone();
            limits.max_retained_bytes = limits.max_retained_bytes.min(
                r.budgets
                    .max_retained_bytes
                    .checked_sub(retained)
                    .filter(|n| *n > 0)
                    .ok_or_else(|| bad("observed_sam_annotation_memory_budget"))?,
            );
            let input = kernel::SamIntakeRequest {
                sam: r.sam.clone(),
                manifest: r.manifest.clone(),
                dictionary: r.dictionary.clone(),
                conversion: r.conversion.clone(),
                reference_id: r.reference_id.clone(),
                reference: r.reference.clone(),
                reference_gene: r.reference_gene.clone(),
                comparison_genes: r.comparison_genes.clone(),
                annotation_sources: sources,
                version_policy: r.version_policy,
                unpadded_reference: r.unpadded_reference,
                required_flank_bases: r.required_flank_bases,
                limits,
                max_units: r.budgets.max_units,
                kernel_budgets: r.kernel_budgets.clone(),
            };
            let (analysis, intake) = kernel::summarize_sam_bound(
                &input,
                &states,
                r.budgets.max_projection_work,
                &mut || ctl.check(),
            )?;
            data.intake = Some(intake);
            data.analysis = Some(analysis);
            ctl.completed = 3;
            if ctl.stop(true) {
                return Err(bad(INTERRUPTED));
            }
            Ok(())
        })();
        if let Err(e) = result {
            if e.message != INTERRUPTED {
                return Err(e);
            }
        }
        loop {
            if ctl.stopped {
                data.analysis = None;
            }
            let stage = staging(
                out.to_str()
                    .ok_or_else(|| bad("observed_sam_output_path"))?,
            )?;
            let status = ctl.status().to_owned();
            let artifacts = publication::artifacts(stage.path(), &data, &status, &r, &mut || {
                if ctl.stopped { Ok(()) } else { ctl.check() }
            });
            let artifacts = match artifacts {
                Ok(a) => a,
                Err(e) if e.message == INTERRUPTED && ctl.stopped => {
                    drop(stage);
                    continue;
                }
                Err(e) => return Err(e),
            };
            ctl.completed = 4;
            if !ctl.stopped && ctl.stop(true) {
                drop(stage);
                continue;
            }
            let status = match ctl.status() {
                "completed" => ObservedSamExecutionStatus::Completed,
                "cancelled" => ObservedSamExecutionStatus::Cancelled,
                _ => ObservedSamExecutionStatus::ControlError,
            };
            let mut report = ObservedSamReport {
                schema: OBSERVED_SAM_REPORT_SCHEMA.into(),
                content_sha256: String::new(),
                request_sha256: identity.clone(),
                policy_sha256: policy(&r)?,
                semantic_context_sha256: data.analysis.as_ref().map(|a| a.context_sha256.clone()),
                request: r.clone(),
                resolved_output_dir: out.display().to_string(),
                status,
                control_error: ctl.error.clone(),
                source_eof_verified: data.intake.as_ref().is_some_and(|i| i.source_eof_verified),
                coordinate_authority: ObservedSamCoordinateAuthority::Declared,
                reference_verification: None,
                identity_scope: "sam_declared_query_unit".into(),
                counts: data.analysis.as_ref().map(report_counts),
                gene_resolutions: data
                    .analysis
                    .as_ref()
                    .map_or_else(Vec::new, |a| a.gene_resolutions.clone()),
                artifacts,
                nonclaims: nonclaims(),
            };
            size(&report, (r.budgets.max_retained_bytes as u64 - 65536) / 64)?;
            report.content_sha256 = digest(&report)?;
            validate_report(&report)?;
            let publication = publication::publish(stage, &report, &mut |_phase| {
                if !ctl.stopped && ctl.stop(true) {
                    Err(bad(INTERRUPTED))
                } else {
                    Ok(())
                }
            });
            match publication {
                Err(e) if e.message == INTERRUPTED && ctl.stopped => continue,
                Err(e) => return Err(e),
                Ok(()) => {
                    let _ = (ctl.progress)(OperationProgress::Workflow {
                        completed: 5,
                        total: 5,
                    });
                    return Ok(report);
                }
            }
        }
    }

    pub fn observed_sam_result_binding(
        request: &Request,
        report: &ObservedSamReport,
    ) -> Result<ObservedSamResultBinding, EngineError> {
        if digest(request)? != report.request_sha256
            || absolute(&request.output_dir)?.to_string_lossy() != report.resolved_output_dir
        {
            return Err(bad("observed_sam_executed_request_mismatch"));
        }
        validate_report(report)
    }
}

fn validate_report(r: &ObservedSamReport) -> Result<ObservedSamResultBinding, EngineError> {
    admit(&r.request)?;
    size(
        r,
        (r.request.budgets.max_retained_bytes as u64 - 65536) / 64,
    )?;
    let mut basis = r.clone();
    basis.content_sha256.clear();
    let names = r
        .artifacts
        .iter()
        .map(|a| a.path.as_str())
        .collect::<BTreeSet<_>>();
    let expected = TABLES.into_iter().chain([DATA]).collect::<BTreeSet<_>>();
    let completed = r.status == ObservedSamExecutionStatus::Completed;
    let mut resolution_ids = BTreeSet::new();
    let resolutions_valid = r.gene_resolutions.iter().all(|row| {
        let Some(source) = r
            .request
            .annotation_sources
            .iter()
            .find(|s| s.selector().source_id == row.source_id)
        else {
            return false;
        };
        resolution_ids.insert(&row.source_id)
            && match row.availability {
                ObservedSamAnnotationAvailability::Available => {
                    row.comparison.as_ref().is_some_and(|c| {
                        c.source_id == row.source_id
                            && c.gene_id.as_deref() == Some(source.selector().gene_id.as_str())
                    })
                }
                ObservedSamAnnotationAvailability::MissingGene
                | ObservedSamAnnotationAvailability::Refused => row.comparison.is_none(),
            }
    });
    if r.schema != OBSERVED_SAM_REPORT_SCHEMA
        || digest(&basis)? != r.content_sha256
        || digest(&r.request)? != r.request_sha256
        || r.policy_sha256 != policy(&r.request)?
        || r.identity_scope != "sam_declared_query_unit"
        || r.coordinate_authority != ObservedSamCoordinateAuthority::Declared
        || r.reference_verification.is_some()
        || !Path::new(&r.resolved_output_dir).is_absolute()
        || kernel::validate_sam_raw_path(&r.resolved_output_dir).is_err()
        || (r.status == ObservedSamExecutionStatus::ControlError) != r.control_error.is_some()
        || r.control_error.as_ref().is_some_and(|s| s.is_empty())
        || completed != r.counts.is_some()
        || completed != r.semantic_context_sha256.is_some()
        || (completed && !r.source_eof_verified)
        || r.semantic_context_sha256
            .as_ref()
            .is_some_and(|s| !is_hash(s))
        || (completed && r.gene_resolutions.len() != r.request.annotation_sources.len())
        || (completed && !resolutions_valid)
        || (!completed && !r.gene_resolutions.is_empty())
        || names != expected
        || r.artifacts.len() != expected.len()
        || r.artifacts
            .iter()
            .try_fold(0u64, |sum, a| {
                if is_hash(&a.sha256) {
                    sum.checked_add(a.bytes)
                } else {
                    None
                }
            })
            .is_none_or(|n| n > r.request.budgets.max_output_bytes)
        || r.counts.as_ref().is_some_and(|c| {
            c.assessed_units.checked_add(c.unavailable_units) != Some(c.query_units)
                || c.query_units > r.request.budgets.max_units
                || c.observations > r.request.text_budgets.max_rows
                || c.optimal_observations > c.observations
                || c.chain_families > r.request.kernel_budgets.max_chain_families
        })
    {
        return Err(bad("observed_sam_report_binding"));
    }
    Ok(ObservedSamResultBinding {
        schema: OBSERVED_SAM_BINDING_SCHEMA.into(),
        request_sha256: r.request_sha256.clone(),
        report_content_sha256: r.content_sha256.clone(),
        report_path: Path::new(&r.resolved_output_dir)
            .join("report.json")
            .display()
            .to_string(),
        outcome: match r.status {
            ObservedSamExecutionStatus::Completed => "completed",
            ObservedSamExecutionStatus::Cancelled => "cancelled",
            ObservedSamExecutionStatus::ControlError => "control_error",
        }
        .into(),
    })
}
