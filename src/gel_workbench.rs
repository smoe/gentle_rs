//! Shared measured-gel workbench. Numerical distributions are conditional on
//! explicit common coverage and their declared integration measure, not counts
//! of independent molecules. Signed residuals are not probability measures.

pub mod import;
#[cfg(test)]
mod tests;

use crate::{digest_utils::sha256_hex_bytes, gel_image::validate_label};
use base64::{Engine as _, engine::general_purpose::STANDARD};
use gentle_protocol::{EngineError, gel_image::*, gel_workbench::*};

pub const MAX_PROFILE_BINS: usize = 20_000;

pub fn import_profile_csv(
    request: &GelProfileFileRequest,
    bytes: &[u8],
) -> Result<(GelProfile, GelNumericSource), EngineError> {
    if bytes.is_empty() || bytes.len() > 4 * 1024 * 1024 {
        return Err(EngineError::invalid_input(
            "Integrated-bin CSV must be nonempty and at most 4 MiB",
        ));
    }
    let mut csv = csv::ReaderBuilder::new().flexible(false).from_reader(bytes);
    let headers = csv
        .headers()
        .map_err(|e| EngineError::invalid_input(e.to_string()))?;
    if headers.iter().collect::<Vec<_>>() != vec!["left_edge", "right_edge", "weight"] {
        return Err(EngineError::invalid_input(
            "CSV requires exactly left_edge,right_edge,weight (integrated bin weights, not arbitrary point-sampled FU values)",
        ));
    }
    let mut edges = Vec::new();
    let mut weights = Vec::new();
    for row in csv.records() {
        let row = row.map_err(|e| EngineError::invalid_input(e.to_string()))?;
        let number = |s: &str| {
            s.parse::<f64>().map_err(|_| {
                EngineError::invalid_input("CSV contains a nonnumeric bin edge/weight")
            })
        };
        let a = number(&row[0])?;
        let b = number(&row[1])?;
        if !a.is_finite()
            || !b.is_finite()
            || a <= 0.0
            || a >= b
            || edges.last().is_some_and(|last| *last > a)
        {
            return Err(EngineError::invalid_input(
                "CSV bins must be positive, finite, ordered and nonoverlapping",
            ));
        }
        if edges.is_empty() {
            edges.push(a);
        } else if edges.last() != Some(&a) {
            edges.push(a);
            weights.push(None);
        }
        edges.push(b);
        weights.push(if row[2].trim().is_empty() {
            None
        } else {
            Some(number(&row[2])?)
        });
        if weights.len() > MAX_PROFILE_BINS {
            return Err(EngineError::invalid_input(
                "CSV exceeds 20000 bins including explicit gaps",
            ));
        }
    }
    let digest = sha256_hex_bytes(bytes);
    let p=GelProfile{id:request.profile_id.clone(),label:request.label.clone(),run_id:request.run_id.clone(),lane_id:request.lane_id.clone(),
        source_sha256:digest.clone(),evidence_kind:GelEvidenceKind::VendorTrace,size_kind:request.size_kind,
        signal_units:request.signal_units.clone(),measure:request.measure.clone(),edges,weights,
        processing:vec!["integrated_bin_csv_v1; exact original retained; blank weights and absent intervals are unavailable".into()],
        warnings:vec!["Units, integration measure, calibration and marker masking are explicitly caller-declared; this importer does not convert point samples to bin areas or independently verify the vendor's quantitative calibration.".into()]};
    validate_profile(&p)?;
    let source = GelNumericSource {
        source_name: std::path::Path::new(&request.path)
            .file_name()
            .and_then(|x| x.to_str())
            .unwrap_or("profile.csv")
            .into(),
        sha256: digest,
        original_base64: STANDARD.encode(bytes),
        request: GelProfileFileRequest {
            path: std::path::Path::new(&request.path)
                .file_name()
                .and_then(|x| x.to_str())
                .unwrap_or("profile.csv")
                .into(),
            ..request.clone()
        },
    };
    Ok((p, source))
}

pub fn validate_profile(profile: &GelProfile) -> Result<(), EngineError> {
    for (value, name) in [
        (&profile.id, "profile id"),
        (&profile.label, "profile label"),
        (&profile.run_id, "run id"),
        (&profile.lane_id, "lane id"),
        (&profile.signal_units, "signal units"),
        (&profile.measure, "integration measure"),
    ] {
        validate_label(value, name)?;
    }
    if profile.source_sha256.len() != 64
        || !profile.source_sha256.bytes().all(|b| b.is_ascii_hexdigit())
    {
        return Err(EngineError::invalid_input(
            "A profile requires its source SHA-256",
        ));
    }
    if profile.weights.is_empty()
        || profile.weights.len() > MAX_PROFILE_BINS
        || profile.edges.len() != profile.weights.len() + 1
        || profile.edges.iter().any(|x| !x.is_finite() || *x <= 0.0)
        || profile.edges.windows(2).any(|p| p[0] >= p[1])
        || profile.weights.iter().flatten().any(|x| {
            !x.is_finite() || (profile.evidence_kind != GelEvidenceKind::SignedResidual && *x < 0.0)
        })
    {
        return Err(EngineError::invalid_input(
            "Profile needs 1..20000 ordered positive size bins, finite weights and explicit null gaps; non-residual weights must be nonnegative",
        ));
    }
    if profile.processing.len() > 128
        || profile.warnings.len() > 128
        || profile
            .processing
            .iter()
            .chain(&profile.warnings)
            .any(|s| s.len() > 2048)
    {
        return Err(EngineError::invalid_input(
            "Profile provenance exceeds its budget",
        ));
    }
    Ok(())
}

/// Decode fixed-weight display luminance on the original pixel grid.
/// Pixel signals remain display diagnostics; no detector or concentration units
/// can be inferred from an exported image's contrast or colour lookup table.
fn original_pixels(
    image: &GelImageRecord,
) -> Result<image::ImageBuffer<image::Luma<f32>, Vec<f32>>, EngineError> {
    crate::gel_image::validate_gel_image_record(image)?;
    let bytes = STANDARD
        .decode(&image.original_base64)
        .map_err(|_| EngineError::invalid_input("Invalid image bytes"))?;
    if u64::from(image.descriptor.width) * u64::from(image.descriptor.height)
        > crate::gel_image::MAX_GEL_IMAGE_PIXELS
    {
        return Err(EngineError::invalid_input("Image exceeds 32 megapixels"));
    }
    let mut reader = image::ImageReader::new(std::io::Cursor::new(&bytes))
        .with_guessed_format()
        .map_err(|e| EngineError::invalid_input(e.to_string()))?;
    let mut limits = image::Limits::default();
    limits.max_alloc = Some(256 * 1024 * 1024);
    limits.max_image_width = Some(32768);
    limits.max_image_height = Some(32768);
    reader.limits(limits);
    let decoded = reader
        .decode()
        .map_err(|e| EngineError::invalid_input(e.to_string()))?;
    let mut gray = image::ImageBuffer::new(decoded.width(), decoded.height());
    // Work on stored display channels, not image's color-space conversion: that
    // can change the transfer function and overshoot at white. Use f64 for the
    // fixed Rec.709 weights and native alpha without an RGBA-f32 image copy.
    macro_rules! luminance {
        ($pixels:expr,$scale:expr) => {
            for (gray, pixel) in gray.pixels_mut().zip($pixels.pixels()) {
                let channels = pixel.0.map(|v| f64::from(v) / $scale);
                if channels.iter().any(|v| !v.is_finite() || !(0.0..=1.0).contains(v)) {
                    return Err(EngineError::invalid_input(
                        "Image display channels must be finite and normalized to 0..1; import calibrated numeric bin weights for other signal domains",
                    ));
                }
                let rgb = |r, g, b| (2126.0 * r + 7152.0 * g + 722.0 * b) / 10000.0;
                let value = match channels.as_slice() {
                    [l] => *l,
                    [l, a] => l * a + (1.0 - a),
                    [r, g, b] => rgb(r, g, b),
                    [r, g, b, a] => rgb(r, g, b) * a + (1.0 - a),
                    _ => unreachable!("supported native image channels"),
                };
                *gray = image::Luma([value as f32]);
            }
        };
    }
    match &decoded {
        image::DynamicImage::ImageLuma8(p) => luminance!(p, 255.0),
        image::DynamicImage::ImageLuma16(p) => luminance!(p, 65535.0),
        image::DynamicImage::ImageLumaA8(p) => luminance!(p, 255.0),
        image::DynamicImage::ImageLumaA16(p) => luminance!(p, 65535.0),
        image::DynamicImage::ImageRgb8(p) => luminance!(p, 255.0),
        image::DynamicImage::ImageRgb16(p) => luminance!(p, 65535.0),
        image::DynamicImage::ImageRgba8(p) => luminance!(p, 255.0),
        image::DynamicImage::ImageRgba16(p) => luminance!(p, 65535.0),
        image::DynamicImage::ImageRgb32F(p) => luminance!(p, 1.0),
        image::DynamicImage::ImageRgba32F(p) => luminance!(p, 1.0),
        _ => {
            return Err(EngineError::invalid_input(
                "Unsupported display pixel format",
            ));
        }
    }
    Ok(gray)
}

fn lane_trace(
    pixels: &image::ImageBuffer<image::Luma<f32>, Vec<f32>>,
    lane: &GelImageLane,
    direction: GelMigrationDirection,
    dark: bool,
) -> Result<Vec<(GelImagePoint, f64)>, EngineError> {
    if ![lane.min.x, lane.min.y, lane.max.x, lane.max.y]
        .iter()
        .all(|x| x.is_finite())
        || lane.min.x < 0.0
        || lane.min.y < 0.0
        || lane.max.x < lane.min.x
        || lane.max.y < lane.min.y
        || lane.max.x > f64::from(pixels.width() - 1)
        || lane.max.y > f64::from(pixels.height() - 1)
    {
        return Err(EngineError::invalid_input(
            "Lane must lie inside original image pixel centers",
        ));
    }
    let horizontal = matches!(
        direction,
        GelMigrationDirection::Left | GelMigrationDirection::Right
    );
    let (start, end, across_start, across_end) = if horizontal {
        (
            lane.min.x.ceil() as u32,
            lane.max.x.floor() as u32,
            lane.min.y.ceil() as u32,
            lane.max.y.floor() as u32,
        )
    } else {
        (
            lane.min.y.ceil() as u32,
            lane.max.y.floor() as u32,
            lane.min.x.ceil() as u32,
            lane.max.x.floor() as u32,
        )
    };
    if start >= end || across_start > across_end || end - start > MAX_PROFILE_BINS as u32 {
        return Err(EngineError::invalid_input(
            "Lane needs 2..20001 migration pixels and at least one transverse pixel",
        ));
    }
    let mut trace = Vec::new();
    for along in start..=end {
        let mut sum = 0.0;
        for across in across_start..=across_end {
            let (x, y) = if horizontal {
                (along, across)
            } else {
                (across, along)
            };
            let value = f64::from(pixels.get_pixel(x, y).0[0]);
            sum += if dark { 1.0 - value } else { value };
        }
        let center = if horizontal {
            GelImagePoint {
                x: f64::from(along),
                y: (lane.min.y + lane.max.y) / 2.0,
            }
        } else {
            GelImagePoint {
                x: (lane.min.x + lane.max.x) / 2.0,
                y: f64::from(along),
            }
        };
        trace.push((center, sum / f64::from(across_end - across_start + 1)));
    }
    if matches!(
        direction,
        GelMigrationDirection::Up | GelMigrationDirection::Left
    ) {
        trace.reverse();
    }
    Ok(trace)
}

pub fn detect_bands(
    image: &GelImageRecord,
    request: &GelDetectRequest,
) -> Result<GelBandSuggestions, EngineError> {
    if request.image_id != image.descriptor.image_id
        || request.lanes.is_empty()
        || request.lanes.len() > 256
        || !(0.0 < request.minimum_prominence && request.minimum_prominence <= 1.0)
        || request.minimum_distance_px == 0
        || request.minimum_distance_px > 1000
        || request.smoothing_radius_px > 32
    {
        return Err(EngineError::invalid_input(
            "Band detection requires bound image, 1..256 lanes, prominence (0,1], distance 1..1000 and smoothing radius 0..32",
        ));
    }
    let area = request.lanes.iter().fold(0u64, |total, lane| {
        let width = (lane.max.x - lane.min.x + 1.0).ceil().max(0.0) as u64;
        let height = (lane.max.y - lane.min.y + 1.0).ceil().max(0.0) as u64;
        total.saturating_add(width.saturating_mul(height))
    });
    if area > 64 * 1024 * 1024 {
        return Err(EngineError::invalid_input(
            "Band detection exceeds the 64-megapixel total lane-scan budget",
        ));
    }
    let mut bands = Vec::new();
    let pixels = original_pixels(image)?;
    let mut ids = std::collections::BTreeSet::new();
    for lane in &request.lanes {
        validate_label(&lane.id, "lane id")?;
        if !ids.insert(&lane.id) {
            return Err(EngineError::invalid_input("Duplicate lane id"));
        }
        let trace = lane_trace(&pixels, lane, request.migration, request.dark_bands)?;
        let radius = request.smoothing_radius_px as usize;
        let smooth: Vec<f64> = (0..trace.len())
            .map(|i| {
                let a = i.saturating_sub(radius);
                let b = (i + radius + 1).min(trace.len());
                trace[a..b].iter().map(|p| p.1).sum::<f64>() / (b - a) as f64
            })
            .collect();
        let mut candidates = Vec::new();
        let window = (request.minimum_distance_px as usize * 2).max(radius * 2 + 1);
        let mut i = 1;
        while i + 1 < smooth.len() {
            let start = i;
            while i + 1 < smooth.len() && smooth[i + 1] == smooth[start] {
                i += 1;
            }
            if i + 1 < smooth.len()
                && smooth[start] > smooth[start - 1]
                && smooth[i] > smooth[i + 1]
            {
                let left = smooth[start.saturating_sub(window)..start]
                    .iter()
                    .copied()
                    .fold(smooth[start], f64::min);
                let right = smooth[i + 1..(i + window + 1).min(smooth.len())]
                    .iter()
                    .copied()
                    .fold(smooth[i], f64::min);
                let prominence = smooth[start] - left.max(right);
                if prominence >= request.minimum_prominence {
                    candidates.push(((start + i) / 2, prominence));
                }
            }
            i += 1;
        }
        candidates.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
        let mut accepted: Vec<usize> = Vec::new();
        for (i, _) in candidates {
            if accepted
                .iter()
                .all(|j| i.abs_diff(*j) >= request.minimum_distance_px as usize)
            {
                accepted.push(i);
            }
        }
        accepted.sort_unstable();
        for (n, i) in accepted.into_iter().enumerate() {
            if bands.len() >= 2048 {
                return Err(EngineError::invalid_input(
                    "More than 2048 band suggestions; restrict lanes or raise prominence",
                ));
            }
            bands.push(GelSampleBand {
                id: format!("suggested_{}_{}", lane.id, n + 1),
                lane_id: lane.id.clone(),
                label: format!("Suggested band {} (unconfirmed)", n + 1),
                center: trace[i].0,
                position_half_width_px: None,
            });
        }
    }
    Ok(GelBandSuggestions{image_sha256:image.descriptor.sha256.clone(),request:request.clone(),bands,
        warnings:vec!["Automatic image peaks are unconfirmed proposals, not ladder identities, molecule identities or concentration measurements. Review plateaus, saturation, labels and smears before accepting. Stored channels use fixed Rec.709 display weights without color-space/transfer conversion (stored_rgb_rec709_v1), with native alpha composited over white (alpha_white_v1); source bytes are unchanged.".into()]})
}

pub fn extract_profile(
    image: &GelImageRecord,
    report: &GelImageAnalysisReport,
    request: &GelExtractProfileRequest,
) -> Result<GelProfile, EngineError> {
    if !(0.0..=1.0).contains(&request.background) || report.request.report_id != request.analysis_id
    {
        return Err(EngineError::invalid_input(
            "Background must be finite 0..1 and analysis must match",
        ));
    }
    if request.masked_ranges.len() > 128
        || request
            .masked_ranges
            .iter()
            .any(|[a, b]| !a.is_finite() || !b.is_finite() || *a <= 0.0 || a >= b)
    {
        return Err(EngineError::invalid_input(
            "Masked size ranges require at most128 finite positive increasing intervals",
        ));
    }
    let mut checked = crate::gel_image::analyze_gel_image(&image.descriptor, &report.request)?;
    checked.gentle_version = report.gentle_version.clone();
    if &checked != report {
        return Err(EngineError::invalid_input(
            "Stored calibration was modified",
        ));
    }
    let lane = report
        .request
        .lanes
        .iter()
        .find(|l| l.id == request.lane_id)
        .ok_or_else(|| EngineError::invalid_input("Unknown lane"))?;
    if lane.id == report.request.ladder.lane_id {
        return Err(EngineError::invalid_input(
            "Do not compare the calibration ladder as a sample profile",
        ));
    }
    let trace = lane_trace(
        &original_pixels(image)?,
        lane,
        report.request.migration,
        request.dark_bands,
    )?;
    let mut bins = Vec::new();
    for pair in trace.windows(2) {
        let a = crate::gel_image::interpolate(
            &report.calibration_bands,
            report.request.migration,
            report.request.migration.coordinate(pair[0].0),
        );
        let b = crate::gel_image::interpolate(
            &report.calibration_bands,
            report.request.migration,
            report.request.migration.coordinate(pair[1].0),
        );
        if let (Some((a, _)), Some((b, _))) = (a, b) {
            let weight = ((pair[0].1 - request.background).max(0.0)
                + (pair[1].1 - request.background).max(0.0))
                / 2.0;
            let masked = request
                .masked_ranges
                .iter()
                .any(|[start, end]| b < *end && a > *start);
            bins.push((b, a, (!masked).then_some(weight)));
        }
    }
    bins.reverse();
    if bins.is_empty() {
        return Err(EngineError::invalid_input(
            "Lane has no fully calibrated pixel intervals",
        ));
    }
    let mut edges = vec![bins[0].0];
    edges.extend(bins.iter().map(|b| b.1));
    let profile=GelProfile{id:request.profile_id.clone(),label:lane.label.clone(),run_id:request.run_id.clone(),
        lane_id:lane.id.clone(),source_sha256:image.descriptor.sha256.clone(),evidence_kind:GelEvidenceKind::ImageDisplay,
        size_kind:report.request.size_kind,signal_units:"normalized_display_intensity_pixel".into(),edges,
        weights:bins.iter().map(|b|b.2).collect(),measure:"integrated_original_migration_pixels_v1".into(),
        processing:vec![format!("calibration:{}; image:{}; dark_bands:{}; background:{}; stored_rgb_rec709_v1; alpha_white_v1; uniform redistribution within calibrated size bins; masked_ranges:{:?}",report.request.report_id,image.descriptor.sha256,request.dark_bands,request.background,request.masked_ranges)],
        warnings:vec!["Image-derived display-weighted shape only; contrast, background, saturation and marker exclusion need human review. Not fluorescence, mass, molecule counts or yield.".into()]};
    validate_profile(&profile)?;
    Ok(profile)
}

fn compatible(a: &GelProfile, b: &GelProfile) -> Result<(), EngineError> {
    validate_profile(a)?;
    validate_profile(b)?;
    if a.size_kind != b.size_kind
        || a.signal_units != b.signal_units
        || a.measure != b.measure
        || a.evidence_kind != b.evidence_kind
    {
        return Err(EngineError::invalid_input(
            "Profiles need identical size kind, evidence kind, units and integration measure; image display and vendor traces cannot be pooled",
        ));
    }
    Ok(())
}

fn grid(profiles: &[&GelProfile], range: [f64; 2]) -> Result<Vec<f64>, EngineError> {
    if !range.iter().all(|x| x.is_finite() && *x > 0.0) || range[0] >= range[1] {
        return Err(EngineError::invalid_input(
            "Choose a finite positive increasing size range",
        ));
    }
    let mut edges = vec![range[0], range[1]];
    for p in profiles {
        edges.extend(p.edges.iter().filter(|x| range[0] < **x && **x < range[1]));
    }
    edges.sort_by(f64::total_cmp);
    edges.dedup();
    if edges.len() > MAX_PROFILE_BINS + 1 {
        return Err(EngineError::invalid_input("Common grid exceeds 20000 bins"));
    }
    Ok(edges)
}

fn rebin(p: &GelProfile, edges: &[f64]) -> Vec<Option<f64>> {
    edges
        .windows(2)
        .map(|bin| {
            let i = p.edges.partition_point(|x| *x <= bin[0]).saturating_sub(1);
            if i >= p.weights.len() || bin[0] < p.edges[i] || bin[1] > p.edges[i + 1] {
                return None;
            }
            p.weights[i].map(|w| w * ((bin[1] - bin[0]) / (p.edges[i + 1] - p.edges[i])))
        })
        .collect()
}

fn quantile(edges: &[f64], weights: &[Option<f64>], q: f64, total: f64) -> f64 {
    let target = q * total;
    let mut sum = 0.0;
    for (i, w) in weights.iter().enumerate() {
        if let Some(w) = w {
            if *w > 0.0 && sum + *w >= target {
                return edges[i] + (edges[i + 1] - edges[i]) * ((target - sum) / w).clamp(0.0, 1.0);
            }
            sum += w;
        }
    }
    edges[edges.len() - 1]
}

fn summary(
    edges: &[f64],
    weights: &[Option<f64>],
    threshold: Option<f64>,
) -> Result<GelDistributionSummary, EngineError> {
    let total = weights.iter().flatten().sum::<f64>();
    if !total.is_finite() || total <= 0.0 {
        return Err(EngineError::invalid_input(
            "No positive finite weight on shared valid coverage; no distribution is available",
        ));
    }
    let fraction = threshold.map(|t| {
        weights
            .iter()
            .enumerate()
            .filter_map(|(i, w)| {
                w.map(|w| w * ((t - edges[i]) / (edges[i + 1] - edges[i])).clamp(0.0, 1.0))
            })
            .sum::<f64>()
            / total
    });
    Ok(GelDistributionSummary {
        total_weight: total,
        q10: quantile(edges, weights, 0.1, total),
        median: quantile(edges, weights, 0.5, total),
        q90: quantile(edges, weights, 0.9, total),
        fraction_below_threshold: fraction,
    })
}

pub fn compare(
    left: &GelProfile,
    right: &GelProfile,
    request: &GelCompareRequest,
) -> Result<GelComparisonReport, EngineError> {
    compatible(left, right)?;
    validate_label(&request.comparison_id, "comparison id")?;
    if left.id != request.left_profile
        || right.id != request.right_profile
        || left.evidence_kind == GelEvidenceKind::SignedResidual
    {
        return Err(EngineError::invalid_input(
            "Bind both profiles; signed residuals do not define probability distributions",
        ));
    }
    let edges = grid(&[left, right], request.range)?;
    if request
        .short_threshold
        .is_some_and(|t| !t.is_finite() || t < request.range[0] || t > request.range[1])
    {
        return Err(EngineError::invalid_input(
            "Threshold must be inside the requested size range",
        ));
    }
    let mut a = rebin(left, &edges);
    let mut b = rebin(right, &edges);
    let mut excluded = Vec::new();
    for i in 0..a.len() {
        if a[i].is_none() || b[i].is_none() {
            a[i] = None;
            b[i] = None;
            excluded.push([edges[i], edges[i + 1]]);
        }
    }
    let sa = summary(&edges, &a, request.short_threshold)?;
    let sb = summary(&edges, &b, request.short_threshold)?;
    let (mut ca, mut cb, mut d, mut w1) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    for i in 0..a.len() {
        let previous = ca - cb;
        if let (Some(a), Some(b)) = (a[i], b[i]) {
            ca += a / sa.total_weight;
            cb += b / sb.total_weight;
        }
        let next = ca - cb;
        d = d.max(previous.abs()).max(next.abs());
        let width = edges[i + 1] - edges[i];
        let factor = if previous * next < 0.0 {
            (previous * previous + next * next) / (2.0 * (previous.abs() + next.abs()))
        } else {
            (previous.abs() + next.abs()) / 2.0
        };
        w1 += width * factor;
    }
    if !w1.is_finite() {
        return Err(EngineError::invalid_input(
            "Distance exceeds finite numeric range",
        ));
    }
    let qq = (1..100)
        .map(|i| {
            let p = f64::from(i) / 100.0;
            GelQqPoint {
                probability: p,
                left: quantile(&edges, &a, p, sa.total_weight),
                right: quantile(&edges, &b, p, sb.total_weight),
            }
        })
        .collect();
    let mut warnings=vec!["Descriptive weighted shape comparison, conditional on shared valid coverage; uniform weight within each recorded size bin. KS distance has no p-value: correlated signal bins are not independent fragment observations or biological replicates.".into(),
        "Separate assay/contrast/loading/calibration effects from biological differences. No gene identity, absolute yield, molecule count or additive priming chemistry is established.".into()];
    warnings.extend(left.warnings.clone());
    warnings.extend(right.warnings.clone());
    if !excluded.is_empty() {
        warnings.push("Coverage gaps/markers were excluded from BOTH profiles, not replaced by zero; inspect excluded_intervals before interpreting quantiles.".into());
    }
    Ok(GelComparisonReport {
        schema: "gentle.gel_comparison.v1".into(),
        request: request.clone(),
        left: left.clone(),
        right: right.clone(),
        common_edges: edges,
        left_weights: a,
        right_weights: b,
        excluded_intervals: excluded,
        left_summary: sa,
        right_summary: sb,
        qq,
        ks_distance: d.clamp(0.0, 1.0),
        wasserstein_distance: w1,
        p_value: None,
        warnings,
    })
}

pub fn algebra(
    profiles: &[&GelProfile],
    request: &GelAlgebraRequest,
) -> Result<GelProfile, EngineError> {
    if profiles.is_empty() || profiles.len() > 16 || profiles.len() != request.terms.len() {
        return Err(EngineError::invalid_input(
            "Algebra needs 1..16 explicit profile terms",
        ));
    }
    for (p, t) in profiles.iter().zip(&request.terms) {
        compatible(profiles[0], p)?;
        if p.id != t.profile_id || !t.coefficient.is_finite() || t.coefficient.abs() > 1e6 {
            return Err(EngineError::invalid_input(
                "Invalid algebra operand/coefficient",
            ));
        }
    }
    let edges = grid(profiles, request.range)?;
    let rebinned: Vec<_> = profiles.iter().map(|p| rebin(p, &edges)).collect();
    let weights = (0..edges.len() - 1)
        .map(|i| {
            let mut sum = 0.0;
            for (r, t) in rebinned.iter().zip(&request.terms) {
                if t.coefficient != 0.0 {
                    sum += r[i]? * t.coefficient;
                }
            }
            Some(sum)
        })
        .collect();
    let mut p = profiles[0].clone();
    p.id = request.output_id.clone();
    p.label = format!("Signed algebra {}", request.output_id);
    p.source_sha256 = sha256_hex_bytes(
        &serde_json::to_vec(&(profiles, request))
            .map_err(|e| EngineError::internal(e.to_string()))?,
    );
    p.evidence_kind = GelEvidenceKind::SignedResidual;
    p.edges = edges;
    p.weights = weights;
    p.processing =
        vec![serde_json::to_string(request).map_err(|e| EngineError::internal(e.to_string()))?];
    p.warnings.push("Signed linear combination with explicit coefficients, not a gene-specific enrichment measurement or validated mixture model. Negative residuals retained; no distribution/p-value.".into());
    validate_profile(&p)?;
    Ok(p)
}
