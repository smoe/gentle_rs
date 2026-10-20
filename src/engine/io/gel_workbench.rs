//! Persistent measured-profile orchestration and no-clobber report projection.

use super::*;
use gentle_protocol::gel_workbench::*;
use std::{
    fs::File,
    io::{Read, Write},
    path::Path,
    sync::Arc,
};

impl GentleEngine {
    pub(super) fn gel_retained_source_bytes(&self) -> usize {
        self.state
            .gel_images
            .images
            .values()
            .map(|i| i.descriptor.byte_count)
            .chain(
                self.state
                    .gel_images
                    .workbench
                    .sources
                    .values()
                    .map(|s| s.original_base64.len().saturating_mul(3) / 4),
            )
            .chain(
                self.state
                    .gel_images
                    .workbench
                    .numeric_sources
                    .values()
                    .map(|s| s.original_base64.len().saturating_mul(3) / 4),
            )
            .fold(0usize, usize::saturating_add)
    }

    pub(super) fn apply_gel_workbench(
        &mut self,
        request: GelWorkbenchRequest,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        use crate::gel_workbench as work;
        let mut output = GelWorkbenchResult::default();
        match request {
            GelWorkbenchRequest::ImportSource(request) => {
                if self.state.gel_images.images.contains_key(&request.image_id) {
                    return Err(EngineError::invalid_input(
                        "Choose a new image id; source imports never overwrite",
                    ));
                }
                let mut bytes = Vec::new();
                File::open(&request.path)
                    .and_then(|f| {
                        f.take((crate::gel_image::MAX_GEL_IMAGE_BYTES + 1) as u64)
                            .read_to_end(&mut bytes)
                    })
                    .map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))?;
                let (source, image) = work::import::import_source(&request, &bytes)?;
                let existing = self.gel_retained_source_bytes();
                if existing
                    .saturating_add(bytes.len())
                    .saturating_add(image.descriptor.byte_count)
                    > 128 * 1024 * 1024
                {
                    return Err(EngineError::invalid_input(
                        "Gel sources and images exceed the 128 MiB project budget",
                    ));
                }
                result.gel_image = Some(image.descriptor.clone());
                result.warnings.extend(source.warnings.clone());
                output.source = Some(GelSourceSummary::from(&source));
                self.state
                    .gel_images
                    .images
                    .insert(request.image_id.clone(), Arc::new(image));
                self.state
                    .gel_images
                    .workbench
                    .sources
                    .insert(request.image_id, Arc::new(source));
            }
            GelWorkbenchRequest::DetectBands(request) => {
                let image = self
                    .state
                    .gel_images
                    .images
                    .get(&request.image_id)
                    .ok_or_else(|| {
                        EngineError::new(ErrorCode::NotFound, "Import the gel image first")
                    })?;
                let suggestions = work::detect_bands(image, &request)?;
                result.warnings.extend(suggestions.warnings.clone());
                output.suggestions = Some(suggestions);
            }
            GelWorkbenchRequest::ImportProfile { mut profile } => {
                if profile.evidence_kind == GelEvidenceKind::SignedResidual {
                    return Err(EngineError::invalid_input(
                        "Signed residuals must be created by algebra, not imported; legacy records remain inspectable",
                    ));
                }
                work::validate_profile(&profile)?;
                profile.warnings.push("Imported structured bin weights: source/calibration/measure are caller-declared; original detector calibration is not independently verified by GENtle.".into());
                work::validate_profile(&profile)?;
                self.insert_gel_profile(profile.clone())?;
                output.profile = Some(profile);
            }
            GelWorkbenchRequest::ImportProfileFile(request) => {
                let mut bytes = Vec::new();
                File::open(&request.path)
                    .and_then(|f| f.take(4 * 1024 * 1024 + 1).read_to_end(&mut bytes))
                    .map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))?;
                let (profile, source) = work::import_profile_csv(&request, &bytes)?;
                if self.gel_retained_source_bytes().saturating_add(bytes.len()) > 128 * 1024 * 1024
                {
                    return Err(EngineError::invalid_input(
                        "Gel source retention exceeds 128 MiB",
                    ));
                }
                self.insert_gel_profile(profile.clone())?;
                self.state
                    .gel_images
                    .workbench
                    .numeric_sources
                    .insert(profile.id.clone(), Arc::new(source));
                output.profile = Some(profile);
            }
            GelWorkbenchRequest::ExtractProfile(request) => {
                let report = self.gel_image_analysis(&request.analysis_id)?;
                let image = self
                    .state
                    .gel_images
                    .images
                    .get(&report.image.image_id)
                    .ok_or_else(|| {
                        EngineError::new(ErrorCode::NotFound, "Original image missing")
                    })?;
                let profile = work::extract_profile(image, report, &request)?;
                self.insert_gel_profile(profile.clone())?;
                output.profile = Some(profile);
            }
            GelWorkbenchRequest::Compare(request) => {
                if self
                    .state
                    .gel_images
                    .workbench
                    .comparisons
                    .contains_key(&request.comparison_id)
                    || self.state.gel_images.workbench.comparisons.len() >= 64
                {
                    return Err(EngineError::invalid_input(
                        "Comparison id exists or 64-report project budget reached",
                    ));
                }
                let get = |id: &str| {
                    self.state
                        .gel_images
                        .workbench
                        .profiles
                        .get(id)
                        .ok_or_else(|| {
                            EngineError::new(
                                ErrorCode::NotFound,
                                format!("Profile '{id}' not found"),
                            )
                        })
                };
                let report = work::compare(
                    get(&request.left_profile)?,
                    get(&request.right_profile)?,
                    &request,
                )?;
                result.warnings.extend(report.warnings.clone());
                output.comparison = Some(report.clone());
                self.state
                    .gel_images
                    .workbench
                    .comparisons
                    .insert(request.comparison_id, Arc::new(report));
            }
            GelWorkbenchRequest::Algebra(request) => {
                let profiles = request
                    .terms
                    .iter()
                    .map(|t| {
                        self.state
                            .gel_images
                            .workbench
                            .profiles
                            .get(&t.profile_id)
                            .map(Arc::as_ref)
                            .ok_or_else(|| {
                                EngineError::new(ErrorCode::NotFound, "Algebra profile not found")
                            })
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let profile = work::algebra(&profiles, &request)?;
                self.insert_gel_profile(profile.clone())?;
                output.profile = Some(profile);
            }
            GelWorkbenchRequest::Inspect {} => {}
            GelWorkbenchRequest::InspectProfile { profile_id } => {
                output.profile = Some(
                    self.state
                        .gel_images
                        .workbench
                        .profiles
                        .get(&profile_id)
                        .ok_or_else(|| {
                            EngineError::new(ErrorCode::NotFound, "Saved profile not found")
                        })?
                        .as_ref()
                        .clone(),
                );
            }
            GelWorkbenchRequest::InspectComparison { comparison_id } => {
                output.comparison = Some(
                    self.state
                        .gel_images
                        .workbench
                        .comparisons
                        .get(&comparison_id)
                        .ok_or_else(|| {
                            EngineError::new(ErrorCode::NotFound, "Saved comparison not found")
                        })?
                        .as_ref()
                        .clone(),
                );
            }
            GelWorkbenchRequest::InspectAlgebra { profile_id } => {
                output.algebra = Some(self.gel_algebra_presentation(&profile_id)?);
            }
            GelWorkbenchRequest::Export(request) => {
                let bytes = self.gel_workbench_export_bytes(&request)?;
                let path = Path::new(&request.path);
                let parent = path
                    .parent()
                    .filter(|p| !p.as_os_str().is_empty())
                    .unwrap_or(Path::new("."));
                let mut temp = tempfile::NamedTempFile::new_in(parent)
                    .map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))?;
                temp.write_all(&bytes)
                    .and_then(|_| temp.as_file().sync_all())
                    .map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))?;
                temp.persist_noclobber(path).map_err(|e| {
                    EngineError::new(
                        ErrorCode::Io,
                        format!("Report export never overwrites: {e}"),
                    )
                })?;
                output.exported_sha256 = Some(crate::digest_utils::sha256_hex_bytes(&bytes));
            }
        }
        output.profile_ids = self
            .state
            .gel_images
            .workbench
            .profiles
            .keys()
            .cloned()
            .collect();
        output.comparison_ids = self
            .state
            .gel_images
            .workbench
            .comparisons
            .keys()
            .cloned()
            .collect();
        result.gel_workbench = Some(Box::new(output));
        result.messages.push("Gel workbench request completed; image-derived shape is not molecular abundance or biological acceptance".into());
        Ok(())
    }

    fn insert_gel_profile(&mut self, profile: GelProfile) -> Result<(), EngineError> {
        if self.state.gel_images.workbench.profiles.len() >= 64
            || self
                .state
                .gel_images
                .workbench
                .profiles
                .contains_key(&profile.id)
        {
            return Err(EngineError::invalid_input(
                "Profile id exists or 64-profile project budget reached",
            ));
        }
        self.state
            .gel_images
            .workbench
            .profiles
            .insert(profile.id.clone(), Arc::new(profile));
        Ok(())
    }

    fn gel_algebra_presentation(&self, id: &str) -> Result<GelAlgebraPresentation, EngineError> {
        let mut presentation =
            crate::gel_workbench::inspect_algebra(&self.state.gel_images.workbench, id)?;
        if presentation.verification == GelAlgebraVerification::Verified {
            for profile in &presentation.inputs {
                if profile.evidence_kind == GelEvidenceKind::SignedResidual {
                    continue;
                }
                let verified = (|| -> Result<(), EngineError> {
                    self.verified_gel_numeric_source(profile)?;
                    for image in self
                        .state
                        .gel_images
                        .images
                        .values()
                        .filter(|i| i.descriptor.sha256 == profile.source_sha256)
                    {
                        crate::gel_image::validate_gel_image_record(image)?;
                    }
                    for source in self
                        .state
                        .gel_images
                        .workbench
                        .sources
                        .values()
                        .filter(|s| s.extracted_image_sha256 == profile.source_sha256)
                    {
                        self.verified_gel_source_summary(source)?;
                    }
                    Ok(())
                })();
                if let Err(error) = verified {
                    presentation.verification = GelAlgebraVerification::Mismatch;
                    presentation.reason = Some(error.message);
                    presentation.series = None;
                    presentation.contributions.clear();
                    break;
                }
            }
        }
        Ok(presentation)
    }

    fn verified_gel_source_summary(
        &self,
        source: &GelSourceRecord,
    ) -> Result<GelSourceSummary, EngineError> {
        use base64::{Engine as _, engine::general_purpose::STANDARD};
        if source.original_base64.len() > crate::gel_image::MAX_GEL_IMAGE_BYTES.div_ceil(3) * 4 {
            return Err(EngineError::invalid_input("Oversized retained source"));
        }
        let original = STANDARD
            .decode(&source.original_base64)
            .map_err(|_| EngineError::invalid_input("Invalid stored source encoding"))?;
        let image = self
            .state
            .gel_images
            .images
            .get(&source.image_id)
            .ok_or_else(|| {
                EngineError::new(ErrorCode::NotFound, "Retained source image missing")
            })?;
        crate::gel_image::validate_gel_image_record(image)?;
        if crate::digest_utils::sha256_hex_bytes(&original) != source.source_sha256
            || source.extracted_image_sha256 != image.descriptor.sha256
        {
            return Err(EngineError::invalid_input(
                "Stored source no longer matches its PDF/archive or extracted image digest",
            ));
        }
        let (replayed, rebuilt) =
            crate::gel_workbench::import::import_source(&source.request, &original)?;
        if replayed.image_id != source.image_id
            || replayed.extraction != source.extraction
            || rebuilt.descriptor.sha256 != image.descriptor.sha256
        {
            return Err(EngineError::invalid_input(
                "Retained PDF/archive no longer reproduces its selected raster",
            ));
        }
        Ok(GelSourceSummary::from(source))
    }

    fn verified_gel_numeric_source(
        &self,
        profile: &GelProfile,
    ) -> Result<Option<GelNumericSourceSummary>, EngineError> {
        use base64::{Engine as _, engine::general_purpose::STANDARD};
        let Some(source) = self
            .state
            .gel_images
            .workbench
            .numeric_sources
            .get(&profile.id)
        else {
            return Ok(None);
        };
        if source.original_base64.len() > (4 * 1024 * 1024usize).div_ceil(3) * 4 {
            return Err(EngineError::invalid_input(
                "Oversized retained numeric source",
            ));
        }
        let bytes = STANDARD
            .decode(&source.original_base64)
            .map_err(|_| EngineError::invalid_input("Invalid numeric source encoding"))?;
        if crate::digest_utils::sha256_hex_bytes(&bytes) != source.sha256
            || source.sha256 != profile.source_sha256
        {
            return Err(EngineError::invalid_input(
                "Retained CSV no longer matches its frozen profile source",
            ));
        }
        let (replayed, _) = crate::gel_workbench::import_profile_csv(&source.request, &bytes)?;
        if &replayed != profile {
            return Err(EngineError::invalid_input(
                "Retained CSV and replay metadata no longer reproduce the frozen profile",
            ));
        }
        Ok(Some(GelNumericSourceSummary {
            profile_id: profile.id.clone(),
            source_name: source.source_name.clone(),
            sha256: source.sha256.clone(),
            signal_units: profile.signal_units.clone(),
            measure: profile.measure.clone(),
        }))
    }

    fn gel_workbench_export_bytes(
        &self,
        request: &GelWorkbenchExportRequest,
    ) -> Result<Vec<u8>, EngineError> {
        let extension = Path::new(&request.path)
            .extension()
            .and_then(|x| x.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        if !matches!(extension.as_str(), "pdf" | "json") {
            return Err(EngineError::invalid_input(
                "Workbench export requires a new .pdf or .json file; per-analysis SVG remains available through gel-image export",
            ));
        }
        if request.comparison_display.is_some()
            && (extension != "pdf" || request.comparison_ids.is_empty())
        {
            return Err(EngineError::invalid_input(
                "comparison_display requires PDF export with at least one saved comparison; JSON retains frozen integrals",
            ));
        }
        crate::gel_image::validate_label(&request.title, "report title")?;
        if request.analysis_ids.len() > 8
            || request.comparison_ids.len() > 8
            || request.algebra_ids.len() > 8
            || request.analysis_ids.is_empty()
                && request.comparison_ids.is_empty()
                && request.algebra_ids.is_empty()
        {
            return Err(EngineError::invalid_input(
                "Choose saved analyses, comparisons and/or algebra results (at most eight of each)",
            ));
        }
        let mut analyses = Vec::new();
        let mut previews = Vec::new();
        let mut comparisons = Vec::new();
        let mut sources = Vec::new();
        let mut numeric_sources = Vec::new();
        let mut algebra = Vec::new();
        let mut algebra_sources = Vec::new();
        let mut algebra_numeric_sources = Vec::new();
        for id in &request.algebra_ids {
            let presentation = self.gel_algebra_presentation(id)?;
            if presentation.verification != GelAlgebraVerification::Verified {
                return Err(EngineError::invalid_input(format!(
                    "Algebra export refused ({:?}): {}",
                    presentation.verification,
                    presentation.reason.as_deref().unwrap_or("unverifiable")
                )));
            }
            for profile in &presentation.inputs {
                if let Some(summary) = self.verified_gel_numeric_source(profile)? {
                    if !algebra_numeric_sources.contains(&summary) {
                        algebra_numeric_sources.push(summary);
                    }
                }
                for source in self
                    .state
                    .gel_images
                    .workbench
                    .sources
                    .values()
                    .filter(|s| s.extracted_image_sha256 == profile.source_sha256)
                {
                    let summary = self.verified_gel_source_summary(source)?;
                    if !algebra_sources.contains(&summary) {
                        algebra_sources.push(summary);
                    }
                }
            }
            algebra.push(presentation);
        }
        let algebra_bins: usize = algebra
            .iter()
            .map(|r| {
                r.inputs.iter().map(|p| p.weights.len()).sum::<usize>()
                    + r.series.as_ref().map_or(0, |s| s.weights.len())
                    + r.contributions
                        .iter()
                        .map(|c| c.weights.len())
                        .sum::<usize>()
            })
            .sum();
        if algebra_bins > 1_280_000 {
            return Err(EngineError::invalid_input(
                "Selected algebra exports exceed 1280000 retained input/result bins",
            ));
        }
        for id in &request.analysis_ids {
            let report = self.gel_image_analysis(id)?;
            let image = self
                .state
                .gel_images
                .images
                .get(&report.image.image_id)
                .ok_or_else(|| EngineError::new(ErrorCode::NotFound, "Original image missing"))?;
            let mut verified =
                crate::gel_image::analyze_gel_image(&image.descriptor, &report.request)?;
            verified.gentle_version = report.gentle_version.clone();
            if &verified != report {
                return Err(EngineError::invalid_input(
                    "Stored sizing result is inconsistent",
                ));
            }
            analyses.push(report);
            previews.push(crate::gel_image::gel_image_export_preview(image)?);
            if let Some(source) = self
                .state
                .gel_images
                .workbench
                .sources
                .get(&report.image.image_id)
            {
                let summary = self.verified_gel_source_summary(source)?;
                if !sources.contains(&summary) {
                    sources.push(summary);
                }
            }
        }
        for id in &request.comparison_ids {
            let r = self
                .state
                .gel_images
                .workbench
                .comparisons
                .get(id)
                .ok_or_else(|| EngineError::new(ErrorCode::NotFound, "Saved comparison missing"))?;
            crate::gel_workbench::verify_comparison(r)?;
            for profile in [&r.left, &r.right] {
                for image in self
                    .state
                    .gel_images
                    .images
                    .values()
                    .filter(|i| i.descriptor.sha256 == profile.source_sha256)
                {
                    crate::gel_image::validate_gel_image_record(image)?;
                }
                if let Some(summary) = self.verified_gel_numeric_source(profile)? {
                    if !numeric_sources.contains(&summary) {
                        numeric_sources.push(summary);
                    }
                }
                for source in self
                    .state
                    .gel_images
                    .workbench
                    .sources
                    .values()
                    .filter(|s| s.extracted_image_sha256 == profile.source_sha256)
                {
                    let summary = self.verified_gel_source_summary(source)?;
                    if !sources.contains(&summary) {
                        sources.push(summary);
                    }
                }
            }
            comparisons.push(r.as_ref());
        }
        if extension == "json" {
            for source in &algebra_sources {
                if !sources.contains(source) {
                    sources.push(source.clone());
                }
            }
            for source in &algebra_numeric_sources {
                if !numeric_sources.contains(source) {
                    numeric_sources.push(source.clone());
                }
            }
            let mut report = serde_json::json!({"schema":"gentle.gel_workbench_export.v1","title":request.title,"sources":sources,"numeric_sources":numeric_sources,"analyses":analyses,"comparisons":comparisons});
            if !algebra.is_empty() {
                report["algebra"] = serde_json::to_value(&algebra)
                    .map_err(|e| EngineError::internal(e.to_string()))?;
            }
            return serde_json::to_vec_pretty(&report)
                .map_err(|e| EngineError::internal(e.to_string()));
        }
        let views: Vec<_> = analyses
            .iter()
            .zip(&previews)
            .map(|(a, p)| (*a, p.as_str()))
            .collect();
        let mut pages = if views.is_empty() && comparisons.is_empty() && !algebra.is_empty() {
            vec![]
        } else {
            vec![gentle_render::gel_workbench::render_gel_overview_svg(
                &request.title,
                &views,
            )]
        };
        for chunk in sources.chunks(8) {
            pages.push(gentle_render::gel_workbench::render_sources_svg(chunk));
        }
        for chunk in numeric_sources.chunks(8) {
            pages.push(gentle_render::gel_workbench::render_numeric_sources_svg(
                chunk,
            ));
        }
        for (analysis, preview) in &views {
            pages.push(gentle_render::gel_image::render_gel_image_analysis_svg(
                analysis, preview,
            ));
        }
        for r in comparisons {
            pages.push(
                gentle_render::gel_workbench::render_comparison_svg_with(
                    r,
                    request
                        .comparison_display
                        .unwrap_or(GelComparisonDisplay::BinIntegral),
                )
                .map_err(EngineError::invalid_input)?,
            );
        }
        algebra_sources.retain(|s| !sources.contains(s));
        algebra_numeric_sources.retain(|s| !numeric_sources.contains(s));
        for chunk in algebra_sources.chunks(8) {
            pages.push(gentle_render::gel_workbench::render_sources_svg(chunk));
        }
        for chunk in algebra_numeric_sources.chunks(8) {
            pages.push(gentle_render::gel_workbench::render_numeric_sources_svg(
                chunk,
            ));
        }
        for r in &algebra {
            pages.push(gentle_render::gel_workbench::render_algebra_svg(r));
            pages.extend(
                gentle_render::gel_workbench::render_algebra_provenance_svgs(r)
                    .map_err(EngineError::invalid_input)?,
            );
        }
        if !algebra.is_empty() && pages.len() > 128 {
            return Err(EngineError::invalid_input(
                "Combined algebra PDF exceeds 128 pages; export fewer results or JSON",
            ));
        }
        let temporary =
            tempfile::tempdir().map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))?;
        let mut paths = Vec::new();
        for (i, page) in pages.iter().enumerate() {
            let path = temporary.path().join(format!("page{i}.svg"));
            std::fs::write(&path, page)
                .map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))?;
            paths.push(path);
        }
        let target = temporary.path().join("report.pdf");
        let refs: Vec<_> = paths.iter().map(|p| p.as_path()).collect();
        crate::svg_vector_pdf::render_svg_files_to_vector_pdf(
            &refs,
            &target,
            crate::svg_png::SvgPngRenderOptions::default(),
        )
        .map_err(EngineError::internal)?;
        std::fs::read(target).map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))
    }
}
