//! Local image ingestion and deterministic manual sizing of measured gel bands.
//!
//! No peak finder, image enhancement, intensity quantification, or simulated
//! migration is used here. The caller explicitly confirms ladder assignments;
//! interpolation is confined to their range and works in log10(bp or kDa).

use std::{collections::BTreeSet, io::Cursor, path::Path};

use base64::{Engine as _, engine::general_purpose::STANDARD};
use gentle_protocol::{EngineError, ErrorCode, gel_image::*};
use image::{ImageFormat, ImageReader, Limits};

use crate::digest_utils::sha256_hex_bytes;

#[cfg(test)]
mod tests;

pub const MAX_GEL_IMAGE_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_GEL_IMAGE_PIXELS: u64 = 32 * 1024 * 1024;
const MAX_MARKS: usize = 2048;

pub(crate) fn validate_label(value: &str, name: &str) -> Result<(), EngineError> {
    if value.trim().is_empty() || value.len() > 240 || value.chars().any(char::is_control) {
        return Err(EngineError::invalid_input(format!(
            "{name} must contain 1-240 bytes of non-control text"
        )));
    }
    Ok(())
}

fn image_reader(bytes: &[u8]) -> Result<ImageReader<Cursor<&[u8]>>, EngineError> {
    let mut reader = ImageReader::new(Cursor::new(bytes))
        .with_guessed_format()
        .map_err(|error| {
            EngineError::invalid_input(format!("Cannot identify gel image: {error}"))
        })?;
    if !matches!(
        reader.format(),
        Some(ImageFormat::Png | ImageFormat::Jpeg | ImageFormat::Tiff)
    ) {
        return Err(EngineError::new(
            ErrorCode::Unsupported,
            "Gel images must be PNG, JPEG or TIFF",
        ));
    }
    let mut limits = Limits::default();
    limits.max_image_width = Some(32768);
    limits.max_image_height = Some(32768);
    limits.max_alloc = Some(256 * 1024 * 1024);
    reader.limits(limits);
    Ok(reader)
}

/// Import original bytes without modification. TIFF page zero must be explicitly
/// selected; later pages and automatic multipage/exposure selection are not supported.
pub fn import_gel_image_bytes(
    image_id: &str,
    source_name: &str,
    bytes: &[u8],
    tiff_page: Option<u32>,
) -> Result<GelImageRecord, EngineError> {
    validate_label(image_id, "image_id")?;
    let source_name = Path::new(source_name)
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| EngineError::invalid_input("Image filename must be valid UTF-8"))?;
    validate_label(source_name, "source_name")?;
    if bytes.is_empty() || bytes.len() > MAX_GEL_IMAGE_BYTES {
        return Err(EngineError::invalid_input(
            "Gel image must be nonempty and at most 16 MiB",
        ));
    }
    let reader = image_reader(bytes)?;
    let format = reader.format().expect("checked by image_reader");
    if (format == ImageFormat::Tiff && tiff_page != Some(0))
        || (format != ImageFormat::Tiff && tiff_page.is_some())
    {
        return Err(EngineError::invalid_input(
            "TIFF requires explicit tiff_page=0; only the first page is currently supported. Omit tiff_page for PNG/JPEG.",
        ));
    }
    let (width, height) = image_reader(bytes)?.into_dimensions().map_err(|error| {
        EngineError::invalid_input(format!("Cannot read image dimensions: {error}"))
    })?;
    if width < 2 || height < 2 || u64::from(width) * u64::from(height) > MAX_GEL_IMAGE_PIXELS {
        return Err(EngineError::invalid_input(
            "Gel image must be at least 2x2 pixels and at most 32 megapixels",
        ));
    }
    let decoded = reader
        .decode()
        .map_err(|error| EngineError::invalid_input(format!("Cannot decode gel image: {error}")))?;
    let mut preview = Cursor::new(Vec::new());
    decoded
        .thumbnail(1600, 1600)
        .to_rgb8()
        .write_to(&mut preview, ImageFormat::Png)
        .map_err(|error| EngineError::internal(format!("Cannot create gel preview: {error}")))?;
    let mut warnings = vec![
        "Coordinates use the original decoded pixels (x right, y down); no EXIF rotation or geometric correction is applied.".into(),
        "Preview is display-only; original bytes and bit depth are preserved. Saturation and band quality require user review.".into(),
    ];
    if format == ImageFormat::Jpeg {
        warnings.push(
            "JPEG is lossy; prefer the original PNG or TIFF acquisition for measurements.".into(),
        );
    }
    if format == ImageFormat::Tiff {
        warnings.push("Only explicitly selected TIFF page 0 is analyzed; other pages/exposures are not combined.".into());
    }
    Ok(GelImageRecord {
        descriptor: GelImageDescriptor {
            schema: GEL_IMAGE_SCHEMA.into(),
            image_id: image_id.into(),
            source_name: source_name.into(),
            sha256: sha256_hex_bytes(bytes),
            byte_count: bytes.len(),
            width,
            height,
            format: format!("{format:?}").to_ascii_lowercase(),
            color_type: format!("{:?}", decoded.color()),
            tiff_page,
            warnings,
        },
        original_base64: STANDARD.encode(bytes),
        preview_png_base64: STANDARD.encode(preview.into_inner()),
    })
}

/// Verify persisted original bytes before reusing their geometry or exporting evidence.
pub fn validate_gel_image_record(record: &GelImageRecord) -> Result<(), EngineError> {
    if record.descriptor.schema != GEL_IMAGE_SCHEMA
        || record.original_base64.len() > MAX_GEL_IMAGE_BYTES.div_ceil(3) * 4
    {
        return Err(EngineError::invalid_input(
            "Unsupported or oversized stored gel image",
        ));
    }
    let bytes = STANDARD
        .decode(&record.original_base64)
        .map_err(|_| EngineError::invalid_input("Invalid stored gel-image encoding"))?;
    if bytes.len() != record.descriptor.byte_count
        || sha256_hex_bytes(&bytes) != record.descriptor.sha256
    {
        return Err(EngineError::invalid_input(
            "Stored gel image does not match its source digest/byte count",
        ));
    }
    let dimensions = image_reader(&bytes)?
        .into_dimensions()
        .map_err(|error| EngineError::invalid_input(format!("Invalid stored image: {error}")))?;
    if dimensions != (record.descriptor.width, record.descriptor.height) {
        return Err(EngineError::invalid_input(
            "Stored gel image dimensions do not match the original",
        ));
    }
    Ok(())
}

/// Export images are reconstructed from source bytes, never trusted from an
/// editable project preview. Interactive viewers may cache the imported preview.
pub fn gel_image_export_preview(record: &GelImageRecord) -> Result<String, EngineError> {
    validate_gel_image_record(record)?;
    let bytes = STANDARD
        .decode(&record.original_base64)
        .map_err(|_| EngineError::invalid_input("Invalid stored gel-image encoding"))?;
    let rebuilt = import_gel_image_bytes(
        &record.descriptor.image_id,
        &record.descriptor.source_name,
        &bytes,
        record.descriptor.tiff_page,
    )?;
    Ok(rebuilt.preview_png_base64)
}

fn valid_point(point: GelImagePoint, image: &GelImageDescriptor) -> bool {
    point.x.is_finite()
        && point.y.is_finite()
        && point.x >= 0.0
        && point.y >= 0.0
        && point.x <= f64::from(image.width - 1)
        && point.y <= f64::from(image.height - 1)
}

fn contains(lane: &GelImageLane, point: GelImagePoint) -> bool {
    point.x >= lane.min.x && point.x <= lane.max.x && point.y >= lane.min.y && point.y <= lane.max.y
}

fn interpolate(
    bands: &[GelLadderBand],
    direction: GelMigrationDirection,
    coordinate: f64,
) -> Option<(f64, Vec<String>)> {
    let pair = bands.windows(2).find(|pair| {
        coordinate >= direction.coordinate(pair[0].center)
            && coordinate <= direction.coordinate(pair[1].center)
    })?;
    let a = direction.coordinate(pair[0].center);
    let b = direction.coordinate(pair[1].center);
    let fraction = (coordinate - a) / (b - a);
    let size = 10.0_f64
        .powf(pair[0].size.log10() * (1.0 - fraction) + pair[1].size.log10() * fraction)
        .clamp(pair[1].size, pair[0].size);
    Some((size, pair.iter().map(|band| band.id.clone()).collect()))
}

/// Analyze only caller-confirmed band geometry. This pure path is shared by
/// engine execution, future click-based editors, and numerical regression tests.
pub fn analyze_gel_image(
    image: &GelImageDescriptor,
    request: &GelImageAnalysisRequest,
) -> Result<GelImageAnalysisReport, EngineError> {
    if image.schema != GEL_IMAGE_SCHEMA || image.width < 2 || image.height < 2 {
        return Err(EngineError::invalid_input(
            "Invalid or unsupported gel image descriptor",
        ));
    }
    if request.image_id != image.image_id || request.image_sha256 != image.sha256 {
        return Err(EngineError::invalid_input(
            "Band coordinates are not bound to this image id and SHA-256",
        ));
    }
    validate_label(&request.report_id, "report_id")?;
    validate_label(&request.ladder.label, "ladder label")?;
    validate_label(&request.ladder.source, "ladder source")?;
    if let Some(system) = &request.ladder.gel_system {
        validate_label(system, "gel_system")?;
    }
    if request.lanes.is_empty()
        || request.lanes.len() > 256
        || !(2..=MAX_MARKS).contains(&request.ladder.bands.len())
        || request.sample_bands.is_empty()
        || request.sample_bands.len() > MAX_MARKS
    {
        return Err(EngineError::invalid_input(
            "Require 1-256 lanes, 2-2048 confirmed ladder bands, and 1-2048 sample bands",
        ));
    }
    if request.ladder.prestained
        && (request.size_kind != GelSizeKind::SdsProteinKda || request.ladder.gel_system.is_none())
    {
        return Err(EngineError::invalid_input(
            "Prestained ladders require protein kDa sizing and an explicit gel_system",
        ));
    }
    let mut lane_ids = BTreeSet::new();
    for lane in &request.lanes {
        validate_label(&lane.id, "lane id")?;
        validate_label(&lane.label, "lane label")?;
        if !lane_ids.insert(&lane.id)
            || !valid_point(lane.min, image)
            || !valid_point(lane.max, image)
            || lane.min.x >= lane.max.x
            || lane.min.y >= lane.max.y
        {
            return Err(EngineError::invalid_input(
                "Lane ids must be unique and lane bounds must lie within the original image",
            ));
        }
    }
    let ladder_lane = request
        .lanes
        .iter()
        .find(|lane| lane.id == request.ladder.lane_id)
        .ok_or_else(|| EngineError::invalid_input("Ladder lane does not exist"))?;
    let mut band_ids = BTreeSet::new();
    for band in &request.ladder.bands {
        validate_label(&band.id, "ladder band id")?;
        if !band_ids.insert(&band.id)
            || !valid_point(band.center, image)
            || !contains(ladder_lane, band.center)
            || !band.size.is_finite()
            || band.size <= 0.0
        {
            return Err(EngineError::invalid_input(
                "Ladder bands require unique ids, positive finite sizes, and centers inside their lane",
            ));
        }
    }
    let mut calibration_bands = request.ladder.bands.clone();
    calibration_bands.sort_by(|a, b| {
        request
            .migration
            .coordinate(a.center)
            .total_cmp(&request.migration.coordinate(b.center))
    });
    for pair in calibration_bands.windows(2) {
        if request.migration.coordinate(pair[1].center)
            - request.migration.coordinate(pair[0].center)
            < 1.0
            || pair[0].size <= pair[1].size
        {
            return Err(EngineError::invalid_input(
                "Ladder sizes must strictly decrease along migration; reference centers must be at least one pixel apart",
            ));
        }
    }
    let mut warnings = image.warnings.clone();
    warnings.push("Manual, shared-ladder calibration assumes comparable straight lanes. No smile, tilt, transfer, or perspective correction has been inferred.".into());
    warnings.push("Piecewise interpolation passes through the confirmed references; it is not an independent calibration validation or a concentration measurement.".into());
    if calibration_bands.len() == 2 {
        warnings.push(
            "Only two reference bands: no independent calibration consistency check is possible."
                .into(),
        );
    }
    match request.size_kind {
        GelSizeKind::LinearDnaBp => warnings.push("Size estimates assume linear double-stranded DNA; they do not identify products or size supercoiled/nicked DNA reliably.".into()),
        GelSizeKind::SdsProteinKda => warnings.push("Apparent SDS-PAGE/Western molecular mass only; modifications, gel chemistry and transfer can affect migration. This does not establish protein identity.".into()),
    }
    if request.ladder.prestained {
        warnings.push("Prestained ladder: confirm that the entered apparent sizes are calibrated for the declared gel/buffer system.".into());
    }
    let mut estimates = Vec::new();
    for band in &request.sample_bands {
        validate_label(&band.id, "sample band id")?;
        validate_label(&band.label, "sample band label")?;
        let lane = request
            .lanes
            .iter()
            .find(|lane| lane.id == band.lane_id)
            .ok_or_else(|| EngineError::invalid_input("Sample band refers to a missing lane"))?;
        if !band_ids.insert(&band.id)
            || lane.id == ladder_lane.id
            || !valid_point(band.center, image)
            || !contains(lane, band.center)
        {
            return Err(EngineError::invalid_input(
                "Sample bands require unique ids and centers inside a non-ladder lane",
            ));
        }
        if band.position_half_width_px.is_some_and(|half| {
            !half.is_finite() || half < 0.0 || half > f64::from(image.width.max(image.height))
        }) {
            return Err(EngineError::invalid_input(
                "Band localization half-width must be finite, nonnegative and bounded by the image dimensions",
            ));
        }
        let coordinate = request.migration.coordinate(band.center);
        let interpolation = interpolate(&calibration_bands, request.migration, coordinate);
        let mut band_warnings = vec![];
        let bounds = band.position_half_width_px.and_then(|half| {
            let lo = interpolate(&calibration_bands, request.migration, coordinate - half);
            let hi = interpolate(&calibration_bands, request.migration, coordinate + half);
            match (lo, hi) {
                (Some((larger, _)), Some((smaller, _))) => Some([smaller, larger]),
                _ => {
                    band_warnings.push("Localization span exceeds calibration; no complete localization size bounds reported.".into());
                    None
                }
            }
        });
        if interpolation.is_none() {
            band_warnings.push("Outside calibrated range; no extrapolation performed.".into());
        }
        estimates.push(GelBandSizeEstimate {
            band_id: band.id.clone(),
            lane_id: band.lane_id.clone(),
            label: band.label.clone(),
            center: band.center,
            status: if interpolation.is_some() {
                GelBandSizingStatus::Interpolated
            } else {
                GelBandSizingStatus::OutsideCalibratedRange
            },
            estimated_size: interpolation.as_ref().map(|(size, _)| *size),
            reference_band_ids: interpolation.map(|(_, ids)| ids).unwrap_or_default(),
            localization_size_bounds: bounds,
            warnings: band_warnings,
        });
    }
    Ok(GelImageAnalysisReport {
        schema: GEL_IMAGE_ANALYSIS_SCHEMA.into(),
        algorithm: GEL_CALIBRATION_METHOD.into(),
        gentle_version: env!("CARGO_PKG_VERSION").into(),
        image: image.clone(),
        request: request.clone(),
        calibration_bands,
        estimates,
        warnings,
    })
}

/// A table of measurements with explicit units, reference ids and source digest.
pub fn gel_image_analysis_tsv(report: &GelImageAnalysisReport) -> Result<Vec<u8>, EngineError> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(Vec::new());
    let write_error =
        |error| EngineError::internal(format!("Cannot serialize gel measurements: {error}"));
    writer
        .write_record([
            "report_id",
            "image_sha256",
            "size_kind",
            "unit",
            "lane_id",
            "band_id",
            "label",
            "x_px",
            "y_px",
            "status",
            "estimated_size",
            "reference_band_ids",
            "localization_min",
            "localization_max",
            "warnings",
        ])
        .map_err(write_error)?;
    for band in &report.estimates {
        writer
            .write_record([
                report.request.report_id.clone(),
                report.image.sha256.clone(),
                match report.request.size_kind {
                    GelSizeKind::LinearDnaBp => "linear_dna_bp",
                    GelSizeKind::SdsProteinKda => "sds_protein_kda",
                }
                .into(),
                report.request.size_kind.unit().into(),
                band.lane_id.clone(),
                band.band_id.clone(),
                band.label.clone(),
                band.center.x.to_string(),
                band.center.y.to_string(),
                match band.status {
                    GelBandSizingStatus::Interpolated => "interpolated",
                    GelBandSizingStatus::OutsideCalibratedRange => "outside_calibrated_range",
                }
                .into(),
                band.estimated_size
                    .map(|v| v.to_string())
                    .unwrap_or_default(),
                band.reference_band_ids.join(","),
                band.localization_size_bounds
                    .map(|v| v[0].to_string())
                    .unwrap_or_default(),
                band.localization_size_bounds
                    .map(|v| v[1].to_string())
                    .unwrap_or_default(),
                report
                    .warnings
                    .iter()
                    .chain(band.warnings.iter())
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(" | "),
            ])
            .map_err(write_error)?;
    }
    writer
        .into_inner()
        .map_err(|error| EngineError::internal(error.to_string()))
}
