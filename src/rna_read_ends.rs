//! Authenticated shared read-end service and immutable exports. No adapter-local biology.

use crate::engine::{
    EngineError, GentleEngine, ObservedSamData, ObservedSamExecutionStatus,
    ObservedSamResultBinding, OperationProgress, PrimerPayloadAnnotation, PrimerPayloadFile,
    RnaReadPanelArtifact, SummarizeObservedTranscriptStructuresFromSamRequest,
};
use crate::observed_transcript_structures::{
    self as kernel, indexed_locus as ingress,
    read_ends::{self as core, coverage, evidence, pileup},
};
use serde::{Deserialize, Serialize, de::DeserializeOwned};
use serde_json::{Value, json};
use std::{
    fs::{self, File},
    io::{Read, Write},
    path::{Path, PathBuf},
};

const CAP: u64 = 128 * 1024 * 1024;
const MEMORY: usize = 512 * 1024 * 1024;

#[cfg(test)]
mod tests;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TailInput {
    pub file: PrimerPayloadFile,
    pub raw_signal_manifest: PrimerPayloadFile,
    pub supported_configurations: Vec<evidence::CallerConfiguration>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Sample {
    pub sample_id: String,
    pub analysis: ObservedSamResultBinding,
    pub model: PrimerPayloadFile,
    pub normalized_annotation: PrimerPayloadFile,
    pub composition: evidence::CompositionPolicy,
    pub tail: Option<TailInput>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Build {
    pub schema: String,
    pub source_id: String,
    pub gene_id: String,
    pub transcript_id: String,
    pub samples: Vec<Sample>,
    pub exact_quantile_probabilities: Vec<f64>,
    pub limits: core::Limits,
    pub output_dir: String,
    pub max_output_bytes: u64,
    pub cancel_path: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Format {
    Json,
    Tsv,
    Svg,
    Png,
    Pdf,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "action", rename_all = "snake_case", deny_unknown_fields)]
pub enum Command {
    Build {
        request: Box<Build>,
    },
    Inspect {
        report: PrimerPayloadFile,
        view: pileup::View,
    },
    Export {
        report: PrimerPayloadFile,
        view: pileup::View,
        formats: Vec<Format>,
        output_dir: String,
        max_output_bytes: u64,
        cancel_path: Option<String>,
    },
    SelectLocus {
        request: Box<ingress::Request>,
        analysis: Box<SummarizeObservedTranscriptStructuresFromSamRequest>,
        cancel_path: Option<String>,
    },
    Cancel {
        command: Box<Command>,
    },
}

impl Command {
    pub fn requires_confirmation(&self) -> bool {
        !matches!(self, Self::Inspect { .. })
    }
    pub fn marker(&self) -> Option<&str> {
        match self {
            Self::Build { request } => request.cancel_path.as_deref(),
            Self::Export { cancel_path, .. } | Self::SelectLocus { cancel_path, .. } => {
                cancel_path.as_deref()
            }
            _ => None,
        }
    }
    fn output(&self) -> Option<&str> {
        match self {
            Self::Build { request } => Some(&request.output_dir),
            Self::Export { output_dir, .. } => Some(output_dir),
            Self::SelectLocus { analysis, .. } => Some(&analysis.output_dir),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Completion {
    pub schema: String,
    pub status: String,
    pub request_sha256: String,
    pub content_sha256: String,
    pub output_dir: String,
    pub artifacts: Vec<RnaReadPanelArtifact>,
}

/// Pure executed-request/completion binding for managed late-cancellation retention.
pub(crate) fn completion_binding(
    command: &Command,
    value: &Value,
) -> Result<Completion, EngineError> {
    if !matches!(command, Command::Build { .. } | Command::Export { .. }) {
        return Err(fail("completion_command"));
    }
    let receipt: Completion = serde_json::from_value(value["completion"].clone()).map_err(fail)?;
    let mut basis = receipt.clone();
    basis.content_sha256.clear();
    if receipt.schema != "gentle.rna_read_end_completion.v1"
        || receipt.status != "completed"
        || receipt.request_sha256 != pileup::hash(command)?
        || receipt.content_sha256 != pileup::hash(&basis)?
        || command.output() != Some(receipt.output_dir.as_str())
        || receipt.artifacts.is_empty()
    {
        return Err(fail("completion_request_binding"));
    }
    Ok(receipt)
}

fn fail(reason: impl std::fmt::Display) -> EngineError {
    EngineError::invalid_input(format!("rna_read_end_service_{reason}"))
}
fn json_bytes<T: Serialize>(value: &T) -> Result<Vec<u8>, EngineError> {
    serde_json::to_vec(value).map_err(fail)
}

fn directory(raw: &str) -> Result<PathBuf, EngineError> {
    let p = kernel::validate_sam_raw_path(raw)?;
    if !p.is_absolute() {
        return Err(fail("absolute_output_required"));
    }
    let parent = p.parent().ok_or_else(|| fail("output_parent"))?;
    for ancestor in parent.ancestors() {
        let meta = fs::symlink_metadata(ancestor).map_err(fail)?;
        if !meta.is_dir() || meta.file_type().is_symlink() {
            return Err(fail("output_parent_symlink"));
        }
        #[cfg(windows)]
        {
            use std::os::windows::fs::MetadataExt;
            if meta.file_attributes() & 0x400 != 0 {
                return Err(fail("output_parent_reparse"));
            }
        }
    }
    if fs::symlink_metadata(p).is_ok() {
        return Err(fail("output_exists"));
    }
    Ok(p.into())
}

fn inputs(command: &Command) -> Vec<&PrimerPayloadFile> {
    match command {
        Command::Build { request } => request
            .samples
            .iter()
            .flat_map(|s| {
                std::iter::once(&s.model)
                    .chain(std::iter::once(&s.normalized_annotation))
                    .chain(
                        s.tail
                            .iter()
                            .flat_map(|t| [&t.file, &t.raw_signal_manifest]),
                    )
            })
            .collect(),
        Command::Inspect { report, .. } | Command::Export { report, .. } => vec![report],
        Command::SelectLocus { request, .. } => std::iter::once(&request.alignment)
            .chain([&request.index, &request.dictionary])
            .chain(request.reference_fasta.iter())
            .chain(request.reference_fai.iter())
            .collect(),
        Command::Cancel { .. } => vec![],
    }
}
fn marker_path(command: &Command) -> Result<Option<PathBuf>, EngineError> {
    let Some(raw) = command.marker() else {
        return Ok(None);
    };
    let p = kernel::validate_sam_raw_path(raw)?;
    if !p.is_absolute() {
        return Err(fail("absolute_marker_required"));
    }
    if let Some(out) = command.output() {
        let out = kernel::validate_sam_raw_path(out)?;
        if p == out || p.starts_with(out) {
            return Err(fail("control_collision"));
        }
    }
    // Markers must not replace inputs or result trees, including macOS case aliases.
    let fold = |p: &Path| p.to_string_lossy().to_lowercase();
    for input in inputs(command) {
        let input = kernel::validate_sam_raw_path(&input.path)?;
        if p == input
            || same_file::is_same_file(p, input).unwrap_or(false)
            || ((cfg!(windows) || cfg!(target_os = "macos")) && fold(p) == fold(input))
        {
            return Err(fail("control_collision"));
        }
    }
    if let Command::Build { request } = command {
        for sample in &request.samples {
            let source = Path::new(&sample.analysis.report_path)
                .parent()
                .ok_or_else(|| fail("report_path"))?;
            if p.starts_with(source) {
                return Err(fail("control_collision"));
            }
        }
    }
    Ok(Some(p.into()))
}
fn cancelled(marker: Option<&Path>, identity: &str) -> Result<bool, EngineError> {
    let Some(p) = marker else {
        return Ok(false);
    };
    let m = match fs::symlink_metadata(p) {
        Ok(m) => m,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(false),
        Err(e) => return Err(fail(e)),
    };
    ingress::local_file(p.to_str().ok_or_else(|| fail("marker_path"))?)?;
    if !m.is_file() || m.len() > 128 {
        return Err(fail("marker_invalid"));
    }
    let mut raw = String::new();
    File::open(p)
        .map_err(fail)?
        .take(129)
        .read_to_string(&mut raw)
        .map_err(fail)?;
    if raw.trim() != identity {
        return Err(fail("marker_request_mismatch"));
    }
    Ok(true)
}
fn bound_json<T: DeserializeOwned>(
    file: &PrimerPayloadFile,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<T, EngineError> {
    ingress::local_file(&file.path)?;
    if file.bytes > CAP {
        return Err(fail("input_byte_budget"));
    }
    kernel::read_sam_bound_json(file, MEMORY, check)
}

struct Loaded {
    sample_id: String,
    bindings: coverage::InputBindings,
    analysis: kernel::ObservedSamAnalysis,
    projection: core::Projection,
    model: core::ModelDocument,
    annotation: PrimerPayloadAnnotation,
    composition: evidence::CompositionEvidence,
    tails: Vec<evidence::TailEvidence>,
    tail_file_sha256: Option<String>,
}

fn load(
    sample: &Sample,
    source: &str,
    limits: core::Limits,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Loaded, EngineError> {
    check()?;
    ingress::local_file(&sample.analysis.report_path)?;
    let inspection = GentleEngine::inspect_observed_transcript_structures_from_sam(
        &sample.analysis.report_path,
        Some(&sample.analysis),
        "query_structures.tsv",
        0,
        true,
    )?;
    let report = inspection.report;
    if report.status != ObservedSamExecutionStatus::Completed {
        return Err(fail("analysis_not_completed"));
    }
    let artifact = report
        .artifacts
        .iter()
        .find(|a| a.path == "analysis.json")
        .ok_or_else(|| fail("data_artifact"))?;
    let file = PrimerPayloadFile {
        path: Path::new(&sample.analysis.report_path)
            .parent()
            .ok_or_else(|| fail("report_path"))?
            .join("artifacts/analysis.json")
            .to_string_lossy()
            .into_owned(),
        bytes: artifact.bytes,
        sha256: artifact.sha256.clone(),
    };
    let data: ObservedSamData = bound_json(&file, check)?;
    let analysis = data.analysis.ok_or_else(|| fail("analysis_unavailable"))?;
    let model: core::ModelDocument = bound_json(&sample.model, check)?;
    let annotation: PrimerPayloadAnnotation = bound_json(&sample.normalized_annotation, check)?;
    let binding = analysis
        .annotation_bindings
        .iter()
        .find(|a| a.source_id == source)
        .ok_or_else(|| fail("annotation_source"))?;
    if binding.annotation_sha256 != sample.normalized_annotation.sha256
        || model.annotation_sha256 != sample.normalized_annotation.sha256
        || model.source_id != source
    {
        return Err(fail("annotation_file_binding"));
    }
    let projection = core::project_read_ends(&model, &analysis, &annotation, limits)?;
    check()?;
    let composition = evidence::unavailable_composition(sample.composition.clone())?;
    let tails = if let Some(tail) = &sample.tail {
        ingress::verify(&tail.raw_signal_manifest, CAP, check)?;
        let input: evidence::TailImport = bound_json(&tail.file, check)?;
        evidence::import_tail_estimates(
            &input,
            &tail.supported_configurations,
            &tail.raw_signal_manifest.sha256,
            &analysis,
            limits.max_units,
        )?
    } else {
        vec![]
    };
    let bindings = coverage::InputBindings {
        report_content_sha256: report.content_sha256,
        data_sha256: file.sha256,
        model_sha256: sample.model.sha256.clone(),
        projection_sha256: pileup::hash(&projection)?,
    };
    Ok(Loaded {
        sample_id: sample.sample_id.clone(),
        bindings,
        analysis,
        projection,
        model,
        annotation,
        composition,
        tails,
        tail_file_sha256: sample.tail.as_ref().map(|t| t.file.sha256.clone()),
    })
}

fn completed(
    report: &PrimerPayloadFile,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<pileup::Report, EngineError> {
    let receipt: Completion = bound_json(report, check)?;
    let mut basis = receipt.clone();
    basis.content_sha256.clear();
    if receipt.schema != "gentle.rna_read_end_completion.v1"
        || receipt.status != "completed"
        || receipt.content_sha256 != pileup::hash(&basis)?
    {
        return Err(fail("completion_binding"));
    }
    let path = Path::new(&report.path);
    if path.parent() != Some(Path::new(&receipt.output_dir)) {
        return Err(fail("report_location_binding"));
    }
    let mut data = None;
    for artifact in &receipt.artifacts {
        if artifact.path.contains(['/', '\\']) || artifact.path == "." || artifact.path == ".." {
            return Err(fail("artifact_path"));
        }
        let file = PrimerPayloadFile {
            path: Path::new(&receipt.output_dir)
                .join("artifacts")
                .join(&artifact.path)
                .to_string_lossy()
                .into_owned(),
            sha256: artifact.sha256.clone(),
            bytes: artifact.bytes,
        };
        ingress::verify(&file, CAP, check)?;
        if artifact.path == "pileup.json" {
            data = Some(bound_json(&file, check)?);
        }
    }
    data.ok_or_else(|| fail("pileup_artifact_missing"))
}

fn publish(
    output: &str,
    request_sha256: String,
    contents: Vec<(String, Vec<u8>)>,
    cap: u64,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Completion, EngineError> {
    if !(1..=CAP).contains(&cap) {
        return Err(fail("output_byte_budget"));
    }
    let out = directory(output)?;
    check()?;
    let stage = tempfile::Builder::new()
        .prefix(".gentle-read-ends-")
        .tempdir_in(out.parent().ok_or_else(|| fail("output_parent"))?)
        .map_err(fail)?;
    let root = stage.path().canonicalize().map_err(fail)?;
    let (mut artifacts, mut total) = (vec![], 0u64);
    for (name, bytes) in contents {
        check()?;
        total = total
            .checked_add(bytes.len() as u64)
            .filter(|n| *n <= cap)
            .ok_or_else(|| fail("output_byte_budget"))?;
        if name.contains(['/', '\\']) || name == "." || name == ".." {
            return Err(fail("artifact_path"));
        }
        let path = root.join(&name);
        let mut file = File::options()
            .write(true)
            .create_new(true)
            .open(&path)
            .map_err(fail)?;
        file.write_all(&bytes).map_err(fail)?;
        file.sync_all().map_err(fail)?;
        drop(file);
        let bound = ingress::fingerprint(&path, cap, check)?;
        artifacts.push(RnaReadPanelArtifact {
            path: name,
            bytes: bound.bytes,
            sha256: bound.sha256,
            locus_id: None,
        });
    }
    let mut receipt = Completion {
        schema: "gentle.rna_read_end_completion.v1".into(),
        status: "completed".into(),
        request_sha256,
        content_sha256: String::new(),
        output_dir: out.to_string_lossy().into_owned(),
        artifacts,
    };
    receipt.content_sha256 = pileup::hash(&receipt)?;
    let raw = json_bytes(&receipt)?;
    GentleEngine::publish_read_end_bundle(stage, &out, &raw, &receipt.artifacts, cap, check)?;
    Ok(receipt)
}

/// Shared execution. Hosts must apply request-aware approval before this method touches files.
pub fn execute(
    engine: &GentleEngine,
    command: Command,
    progress: &mut dyn FnMut(OperationProgress) -> bool,
) -> Result<Value, EngineError> {
    if let Command::Cancel { command: target } = &command {
        if matches!(
            target.as_ref(),
            Command::Cancel { .. } | Command::Inspect { .. }
        ) {
            return Err(fail("cancel_target"));
        }
        let identity = pileup::hash(target)?;
        let path = marker_path(target)?.ok_or_else(|| fail("cancel_path_required"))?;
        if cancelled(Some(&path), &identity)? {
            return Ok(json!({"status":"cancellation_requested","termination_confirmed":false}));
        }
        let parent = path.parent().ok_or_else(|| fail("marker_parent"))?;
        let mut staged = tempfile::NamedTempFile::new_in(parent).map_err(fail)?;
        staged.write_all(identity.as_bytes()).map_err(fail)?;
        staged.as_file().sync_all().map_err(fail)?;
        marker_path(target)?;
        staged.persist_noclobber(path).map_err(fail)?;
        return Ok(json!({"status":"cancellation_requested","termination_confirmed":false}));
    }
    let identity = pileup::hash(&command)?;
    let marker = marker_path(&command)?;
    let mut stopped = false;
    let mut check = || {
        stopped |= !progress(OperationProgress::Workflow {
            completed: 0,
            total: 1,
        });
        if stopped || cancelled(marker.as_deref(), &identity)? {
            return Err(fail("cancelled"));
        }
        Ok(())
    };
    check()?;
    match command {
        Command::Build { request } => {
            if request.schema != "gentle.rna_read_end_build.v1"
                || request.samples.is_empty()
                || request.samples.len() > 64
                || request.limits.max_units > 100_000
                || request.limits.max_observations > 1_000_000
                || request.limits.max_work > 100_000_000
            {
                return Err(fail("build_admission"));
            }
            directory(&request.output_dir)?;
            let input_bytes = inputs(&Command::Build {
                request: request.clone(),
            })
            .into_iter()
            .try_fold(0u64, |n, f| n.checked_add(f.bytes))
            .ok_or_else(|| fail("input_byte_budget"))?;
            if input_bytes
                .checked_mul(64)
                .and_then(|n| n.checked_add(65536))
                .is_none_or(|n| n > MEMORY as u64)
            {
                return Err(fail("input_memory_budget"));
            }
            let loaded = request
                .samples
                .iter()
                .map(|s| load(s, &request.source_id, request.limits, &mut check))
                .collect::<Result<Vec<_>, _>>()?;
            let inputs: Vec<_> = loaded
                .iter()
                .map(|l| pileup::Input {
                    sample: coverage::SampleInput {
                        sample_id: &l.sample_id,
                        gene_id: &request.gene_id,
                        bindings: &l.bindings,
                        analysis: &l.analysis,
                        projection: &l.projection,
                    },
                    model: &l.model,
                    annotation: &l.annotation,
                    composition: &l.composition,
                    tails: &l.tails,
                    tail_file_sha256: l.tail_file_sha256.as_deref(),
                })
                .collect();
            let report = pileup::build(
                &request.source_id,
                &request.gene_id,
                &request.transcript_id,
                &request.exact_quantile_probabilities,
                &inputs,
                request.limits,
            )?;
            check()?;
            let contents = vec![
                ("pileup.json".into(), json_bytes(&report)?),
                ("coverage.json".into(), json_bytes(&report.coverage)?),
            ];
            let receipt = publish(
                &request.output_dir,
                identity,
                contents,
                request.max_output_bytes,
                &mut check,
            )?;
            Ok(json!({"completion":receipt,"coverage":report.coverage}))
        }
        Command::Inspect { report, view } => {
            let report = completed(&report, &mut check)?;
            let page = pileup::page(&report, &view)?;
            let scene = pileup::scene(&report, &page)?;
            Ok(json!({"page":page,"scene":scene}))
        }
        Command::Export {
            report,
            view,
            formats,
            output_dir,
            max_output_bytes,
            ..
        } => {
            directory(&output_dir)?;
            if formats.is_empty()
                || formats.len() > 5
                || formats
                    .iter()
                    .enumerate()
                    .any(|(i, f)| formats[..i].contains(f))
            {
                return Err(fail("export_formats"));
            }
            let report = completed(&report, &mut check)?;
            let page = pileup::page(&report, &view)?;
            let scene = pileup::scene(&report, &page)?;
            let svg = gentle_render::read_end_pileup::render(&scene).map_err(fail)?;
            let mut contents = vec![];
            for format in formats {
                check()?;
                let (name, bytes) = match format {
                    Format::Json => ("page.json", json_bytes(&page)?),
                    Format::Tsv => ("lanes.tsv", pileup::tsv(&page)?.into_bytes()),
                    Format::Svg => ("pileup.svg", svg.as_bytes().to_vec()),
                    Format::Png | Format::Pdf => {
                        // Prevent unbounded raster allocation; exact SVG/JSON pages can be larger.
                        let height = 30 + 19 * scene.header.len() + 95 + 54 * scene.rows.len();
                        if (scene.width as usize)
                            .checked_mul(height)
                            .is_none_or(|pixels| pixels > 32_000_000)
                        {
                            return Err(fail("raster_pixel_budget"));
                        }
                        let options = crate::svg_png::SvgPngRenderOptions::default();
                        if format == Format::Png {
                            (
                                "pileup.png",
                                crate::svg_png::render_svg_to_png_bytes(&svg, options)
                                    .map_err(fail)?
                                    .bytes,
                            )
                        } else {
                            (
                                "pileup.pdf",
                                crate::svg_vector_pdf::render_svg_to_vector_pdf_bytes_audited(
                                    &svg, options,
                                )
                                .map_err(fail)?
                                .bytes,
                            )
                        }
                    }
                };
                contents.push((name.into(), bytes));
            }
            let receipt = publish(
                &output_dir,
                identity,
                contents,
                max_output_bytes,
                &mut check,
            )?;
            Ok(json!({"completion":receipt,"view_sha256":page.view_sha256}))
        }
        Command::SelectLocus {
            request,
            mut analysis,
            ..
        } => {
            directory(&analysis.output_dir)?;
            if request.reference_id != analysis.reference_id
                || request.reference.assembly != analysis.reference.assembly
                || request.reference.fasta_sha256 != analysis.reference.fasta_sha256
                || pileup::hash(&request.dictionary)? != pileup::hash(&analysis.dictionary)?
                || analysis.cancel_path.is_some()
            {
                return Err(fail("selection_analysis_binding"));
            }
            let selected = ingress::select(&request, &mut check)?;
            analysis.sam = selected.sam.clone();
            analysis.manifest = selected.manifest.clone();
            analysis.conversion = kernel::ObservedSamConversion::ExecutedIndexedSelection {
                receipt: Box::new(selected.receipt.clone()),
            };
            let report =
                engine.summarize_observed_sam_internal(*analysis, &mut |_| check().is_ok())?;
            if report.status != ObservedSamExecutionStatus::Completed {
                return Err(fail("selection_analysis_not_completed"));
            }
            Ok(
                json!({"observed_sam_report":report,"selection":selected.receipt,"command_sha256":identity}),
            )
        }
        Command::Cancel { .. } => Err(fail("cancel_target")),
    }
}
