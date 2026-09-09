//! Measured-gel operations, immutable project storage, and no-clobber exports.

use std::{
    fs::File,
    io::{Read, Write},
    path::Path,
    sync::Arc,
};

use super::*;
use crate::gel_image::{
    MAX_GEL_IMAGE_BYTES, analyze_gel_image, gel_image_analysis_tsv, import_gel_image_bytes,
    validate_gel_image_record,
};
use gentle_protocol::gel_image::*;

impl GentleEngine {
    /// Retrieve an immutable analysis for GUI/CLI inspection without decoding an image.
    pub fn gel_image_analysis(
        &self,
        report_id: &str,
    ) -> Result<&GelImageAnalysisReport, EngineError> {
        self.state
            .gel_images
            .analyses
            .get(report_id)
            .map(Arc::as_ref)
            .ok_or_else(|| {
                EngineError::new(
                    ErrorCode::NotFound,
                    format!("Gel image analysis '{report_id}' not found"),
                )
            })
    }

    pub(super) fn apply_gel_image_operation(
        &mut self,
        op: Operation,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        match op {
            Operation::ImportGelImage { request } => {
                if self.state.gel_images.images.contains_key(&request.image_id) {
                    return Err(EngineError::invalid_input(
                        "Image id already exists; choose a new id rather than replace the original",
                    ));
                }
                let mut bytes = Vec::new();
                File::open(&request.path)
                    .and_then(|file| {
                        file.take((MAX_GEL_IMAGE_BYTES + 1) as u64)
                            .read_to_end(&mut bytes)
                    })
                    .map_err(|error| {
                        EngineError::new(ErrorCode::Io, format!("Cannot read gel image: {error}"))
                    })?;
                let record = import_gel_image_bytes(
                    &request.image_id,
                    &request.path,
                    &bytes,
                    request.tiff_page,
                )?;
                let retained_bytes = self
                    .state
                    .gel_images
                    .images
                    .values()
                    .map(|image| image.descriptor.byte_count)
                    .fold(0usize, usize::saturating_add);
                if retained_bytes.saturating_add(bytes.len()) > 128 * 1024 * 1024 {
                    return Err(EngineError::invalid_input(
                        "This initial self-contained gel store is limited to 128 MiB of original images per project",
                    ));
                }
                result.warnings.extend(record.descriptor.warnings.clone());
                result.messages.push(format!(
                    "Imported measured gel '{}' ({} x {} pixels); original retained unchanged",
                    request.image_id, record.descriptor.width, record.descriptor.height
                ));
                result.gel_image = Some(record.descriptor.clone());
                self.state
                    .gel_images
                    .images
                    .insert(request.image_id, Arc::new(record));
            }
            Operation::AnalyzeGelImage { request } => {
                if self
                    .state
                    .gel_images
                    .analyses
                    .contains_key(&request.report_id)
                {
                    return Err(EngineError::invalid_input(
                        "Report id already exists; use a new id for revised band assignments",
                    ));
                }
                let image = self
                    .state
                    .gel_images
                    .images
                    .get(&request.image_id)
                    .ok_or_else(|| {
                        EngineError::new(
                            ErrorCode::NotFound,
                            "Import the gel image before analyzing bands",
                        )
                    })?;
                validate_gel_image_record(image)?;
                let report = analyze_gel_image(&image.descriptor, &request)?;
                result.warnings.extend(report.warnings.clone());
                result.messages.push(format!(
                    "Sized {} manually identified bands in {} using confirmed ladder references",
                    report.estimates.len(),
                    request.size_kind.unit()
                ));
                result.gel_image_analysis = Some(Box::new(report.clone()));
                self.state
                    .gel_images
                    .analyses
                    .insert(request.report_id.clone(), Arc::new(report));
            }
            Operation::InspectGelImageAnalysis { report_id } => {
                result.gel_image_analysis =
                    Some(Box::new(self.gel_image_analysis(&report_id)?.clone()));
                result
                    .messages
                    .push(format!("Inspected measured gel analysis '{report_id}'"));
            }
            Operation::ExportGelImageAnalysis { request } => {
                let report = self.gel_image_analysis(&request.report_id)?;
                let image = self
                    .state
                    .gel_images
                    .images
                    .get(&report.image.image_id)
                    .ok_or_else(|| {
                        EngineError::new(ErrorCode::NotFound, "Original gel image is missing")
                    })?;
                validate_gel_image_record(image)?;
                // Revalidate persisted derived values before labelling an exported figure as evidence.
                let mut recomputed = analyze_gel_image(&image.descriptor, &report.request)?;
                recomputed.gentle_version = report.gentle_version.clone();
                if &recomputed != report {
                    return Err(EngineError::invalid_input(
                        "Stored gel analysis is inconsistent with its image/calibration; analyze again with a new report id",
                    ));
                }
                let bytes = match request.format {
                    GelImageExportFormat::Json => serde_json::to_vec_pretty(report)
                        .map_err(|error| EngineError::internal(error.to_string()))?,
                    GelImageExportFormat::Tsv => gel_image_analysis_tsv(report)?,
                    GelImageExportFormat::Svg => {
                        let preview = crate::gel_image::gel_image_export_preview(image)?;
                        gentle_render::gel_image::render_gel_image_analysis_svg(report, &preview)
                            .into_bytes()
                    }
                };
                let path = Path::new(&request.path);
                let parent = path
                    .parent()
                    .filter(|parent| !parent.as_os_str().is_empty())
                    .unwrap_or_else(|| Path::new("."));
                let mut file = tempfile::NamedTempFile::new_in(parent).map_err(|error| {
                    EngineError::new(ErrorCode::Io, format!("Cannot prepare gel export: {error}"))
                })?;
                file.write_all(&bytes)
                    .and_then(|_| file.as_file().sync_all())
                    .map_err(|error| {
                        EngineError::new(ErrorCode::Io, format!("Cannot write gel export: {error}"))
                    })?;
                file.persist_noclobber(path)
                    .map_err(|error| EngineError::new(ErrorCode::Io, format!("Cannot commit gel export (existing files are never overwritten): {error}")))?;
                result.messages.push(format!(
                    "Exported measured gel analysis '{}' to '{}'",
                    request.report_id, request.path
                ));
                result.gel_image_analysis = Some(Box::new(report.clone()));
            }
            _ => unreachable!("not a gel-image operation"),
        }
        Ok(())
    }
}
