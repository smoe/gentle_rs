//! Authenticated locus actions reuse the read-end service's file, control and publication guards.

use super::*;
use crate::observed_transcript_structures::read_ends::locus as kernel;

pub const BUILD_SCHEMA: &str = "gentle.rna_read_end_locus_build.v1";
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Sample {
    pub sample_id: String,
    pub analysis: ObservedSamResultBinding,
    pub annotation: Option<PrimerPayloadFile>,
    pub model: Option<PrimerPayloadFile>,
    pub original_bases: Option<PrimerPayloadFile>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Build {
    pub schema: String,
    pub source_id: String,
    pub region: kernel::Region,
    pub samples: Vec<Sample>,
    pub terminal_runs: oriented::Policy,
    pub limits: core::Limits,
    pub output_dir: String,
    pub max_output_bytes: u64,
    pub cancel_path: Option<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "action", rename_all = "snake_case", deny_unknown_fields)]
pub enum Command {
    Build {
        request: Box<Build>,
    },
    Inspect {
        report: PrimerPayloadFile,
        view: kernel::View,
    },
    Export {
        report: PrimerPayloadFile,
        view: kernel::View,
        formats: Vec<Format>,
        output_dir: String,
        max_output_bytes: u64,
        cancel_path: Option<String>,
    },
}
impl Command {
    pub fn is_inspection(&self) -> bool {
        matches!(self, Self::Inspect { .. })
    }
    pub fn is_build(&self) -> bool {
        matches!(self, Self::Build { .. })
    }
    pub fn marker(&self) -> Option<&str> {
        match self {
            Self::Build { request } => request.cancel_path.as_deref(),
            Self::Export { cancel_path, .. } => cancel_path.as_deref(),
            _ => None,
        }
    }
    pub fn output(&self) -> Option<&str> {
        match self {
            Self::Build { request } => Some(&request.output_dir),
            Self::Export { output_dir, .. } => Some(output_dir),
            _ => None,
        }
    }
    pub(super) fn inputs(&self) -> Vec<&PrimerPayloadFile> {
        match self {
            Self::Build { request } => request
                .samples
                .iter()
                .flat_map(|s| {
                    s.annotation
                        .iter()
                        .chain(s.model.iter())
                        .chain(s.original_bases.iter())
                })
                .collect(),
            Self::Inspect { report, .. } | Self::Export { report, .. } => vec![report],
        }
    }
    pub(super) fn analysis_reports(&self) -> Vec<&str> {
        match self {
            Self::Build { request } => request
                .samples
                .iter()
                .map(|s| s.analysis.report_path.as_str())
                .collect(),
            _ => vec![],
        }
    }
}
struct Loaded {
    sample_id: String,
    bindings: kernel::Bindings,
    analysis: crate::engine::ObservedSamAnalysis,
    annotation: Option<PrimerPayloadAnnotation>,
    model: Option<core::ModelDocument>,
    originals: Option<oriented::Document>,
    original_sha: Option<String>,
    diagnostics: Vec<String>,
}
fn optional<T: DeserializeOwned>(
    file: Option<&PrimerPayloadFile>,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Option<T>, EngineError> {
    let Some(file) = file else {
        return Ok(None);
    };
    crate::observed_transcript_structures::validate_sam_raw_path(&file.path)?;
    match fs::symlink_metadata(&file.path) {
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(e) => Err(fail(e)),
        Ok(_) => bound_json(file, check).map(Some),
    }
}
fn load(
    sample: &Sample,
    remaining: &mut usize,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Loaded, EngineError> {
    check()?;
    let (report, file) = analysis_file(&sample.analysis)?;
    let size = [
        file.bytes,
        fs::metadata(&sample.analysis.report_path)
            .map_err(fail)?
            .len(),
    ]
    .into_iter()
    .chain(
        sample
            .annotation
            .iter()
            .chain(sample.model.iter())
            .chain(sample.original_bases.iter())
            .map(|f| f.bytes),
    )
    .try_fold(0u64, |n, b| n.checked_add(b))
    .and_then(|n| n.checked_mul(64))
    .and_then(|n| n.checked_add(65536))
    .and_then(|n| usize::try_from(n).ok())
    .ok_or_else(|| fail("input_memory_budget"))?;
    *remaining = remaining
        .checked_sub(size)
        .ok_or_else(|| fail("input_memory_budget"))?;
    let data: ObservedSamData = bound_json(&file, check)?;
    let analysis = data.analysis.ok_or_else(|| fail("analysis_unavailable"))?;
    let annotation: Option<PrimerPayloadAnnotation> = optional(sample.annotation.as_ref(), check)?;
    let mut diagnostics = vec![];
    if annotation.is_none() {
        diagnostics.push(if let Some(file) = &sample.annotation {
            format!("annotation_declared_not_found:{}", file.sha256)
        } else {
            "annotation_unavailable_genomic_evidence_retained".into()
        });
    }
    let annotation_sha = annotation
        .as_ref()
        .and(sample.annotation.as_ref())
        .map(|f| f.sha256.clone());
    if annotation_sha.as_ref().is_some_and(|sha| {
        !analysis
            .annotation_bindings
            .iter()
            .any(|b| b.annotation_sha256 == *sha)
    }) {
        diagnostics.push("annotation_changed_scores_recomputed_not_reused".into());
    }
    let supplied_model: Option<core::ModelDocument> = optional(sample.model.as_ref(), check)?;
    if supplied_model.is_none()
        && let Some(file) = &sample.model
    {
        diagnostics.push(format!("cds_model_declared_not_found:{}", file.sha256));
    }
    let model = if annotation.is_some() {
        supplied_model
    } else {
        if supplied_model.is_some() {
            diagnostics.push("cds_without_annotation_not_used".into());
        }
        None
    };
    if model.is_none() && sample.model.is_none() {
        diagnostics.push("cds_models_unavailable_utr_missing".into());
    }
    let originals = sample
        .original_bases
        .as_ref()
        .map(|f| bound_json(f, check))
        .transpose()?;
    Ok(Loaded {
        sample_id: sample.sample_id.clone(),
        bindings: kernel::Bindings {
            report_content_sha256: report.content_sha256,
            data_sha256: file.sha256,
            annotation_sha256: annotation_sha,
            model_sha256: model
                .as_ref()
                .and(sample.model.as_ref())
                .map(|f| f.sha256.clone()),
        },
        analysis,
        annotation,
        model,
        originals,
        original_sha: sample.original_bases.as_ref().map(|f| f.sha256.clone()),
        diagnostics,
    })
}
pub(super) fn execute(
    command: Command,
    identity: String,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Value, EngineError> {
    match command {
        Command::Build { request } => {
            if request.schema != BUILD_SCHEMA
                || request.samples.is_empty()
                || request.samples.len() > 64
            {
                return Err(fail("locus_build_admission"));
            }
            directory(&request.output_dir)?;
            let mut remaining = MEMORY;
            let loaded = request
                .samples
                .iter()
                .map(|s| load(s, &mut remaining, check))
                .collect::<Result<Vec<_>, _>>()?;
            let inputs: Vec<_> = loaded
                .iter()
                .map(|l| kernel::Input {
                    sample_id: &l.sample_id,
                    bindings: l.bindings.clone(),
                    analysis: &l.analysis,
                    annotation: l.annotation.as_ref(),
                    model: l.model.as_ref(),
                    original_bases: l.originals.as_ref().zip(l.original_sha.as_deref()).map(
                        |(document, file_sha256)| oriented::Input {
                            document,
                            file_sha256,
                        },
                    ),
                    diagnostics: l.diagnostics.clone(),
                })
                .collect();
            let report = kernel::build(
                &request.source_id,
                &request.region,
                &inputs,
                request.terminal_runs,
                request.limits,
                check,
            )?;
            let receipt = publish(
                &request.output_dir,
                identity,
                vec![("pileup.json".into(), json_bytes(&report)?)],
                request.max_output_bytes,
                check,
            )?;
            Ok(json!({"report":receipt.file_binding()?,"completion":receipt}))
        }
        Command::Inspect { report, view } => {
            let report: kernel::Report = completed_record(&report, check)?;
            let page = kernel::page(&report, &view)?;
            let scene = kernel::scene(&report, &page)?;
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
            validate_formats(&formats)?;
            let report: kernel::Report = completed_record(&report, check)?;
            let page = kernel::page(&report, &view)?;
            let scene = kernel::scene(&report, &page)?;
            let contents = export_contents(&page, || kernel::tsv(&page), &scene, formats, check)?;
            let receipt = publish(&output_dir, identity, contents, max_output_bytes, check)?;
            Ok(
                json!({"report":receipt.file_binding()?,"completion":receipt,"view_sha256":page.view_sha256}),
            )
        }
    }
}
