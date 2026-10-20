//! Display-only projections of frozen comparison bins. No rebinning, fitting,
//! normalization or replay: callers verify scientific reports before use.

use super::*;

pub(super) fn heading(mode: GelComparisonDisplay) -> &'static str {
    match mode {
        GelComparisonDisplay::BinIntegral => {
            "Raw profile overlay | integrated weight per common bin"
        }
        GelComparisonDisplay::SignalDensityPerSize => {
            "Signal density overlay | saved integral per size width"
        }
    }
}

fn density_bins(report: &GelComparisonReport, bins: &[Bin]) -> Result<Vec<Bin>, &'static str> {
    if report.left.size_kind != report.right.size_kind
        || report.left.signal_units != report.right.signal_units
        || report.left.measure != report.right.measure
        || report.left.evidence_kind != report.right.evidence_kind
        || report.left.evidence_kind == GelEvidenceKind::SignedResidual
    {
        return Err(
            "Density unavailable: incompatible measures/units or signed residual; use algebra inspection for signed data.",
        );
    }
    bins.iter().map(|bin| {
        let mut values = bin.weights;
        for value in values.iter_mut().flatten() {
            let original = *value;
            *value /= bin.end - bin.start;
            if !value.is_finite() || *value < 0.0 || (original > 0.0 && *value == 0.0) {
                return Err("Density unavailable: negative/nonfinite value or overflow; no clipping applied.");
            }
        }
        Ok(Bin { start: bin.start, end: bin.end, weights: values })
    }).collect()
}

pub(super) fn draw_profile_panel(
    svg: &mut Svg,
    report: &GelComparisonReport,
    bins: &[Bin],
    mode: GelComparisonDisplay,
    top: f64,
) -> Result<(), &'static str> {
    let density = mode == GelComparisonDisplay::SignalDensityPerSize;
    let projected;
    let bins = if density {
        projected = density_bins(report, bins)?;
        &projected
    } else {
        bins
    };
    let compatible = report.left.signal_units == report.right.signal_units
        && report.left.measure == report.right.measure
        && report.left.evidence_kind == report.right.evidence_kind;
    let [Some(l), Some(r)] = [weight_range(bins, 0), weight_range(bins, 1)] else {
        return Err("Raw weight range exceeds finite display geometry; overlay omitted.");
    };
    let shared = if density {
        let max = bins
            .iter()
            .flat_map(|b| b.weights)
            .flatten()
            .reduce(f64::max)
            .ok_or("Density unavailable: no evaluated bins.")?;
        [0.0, if max == 0.0 { 1.0 } else { max }]
    } else {
        [l[0].min(r[0]), l[1].max(r[1])]
    };
    let ranges = if compatible { [shared, shared] } else { [l, r] };
    if ranges
        .iter()
        .any(|range| !(range[1] - range[0]).is_finite())
    {
        return Err("Raw weight range exceeds finite display geometry; overlay omitted.");
    }
    let unit = report.left.size_kind.unit();
    if density {
        svg.text(
            32.0,
            top + 22.0,
            "Common-bin weight / bin width; not normalized to unit area; not mass/molar density.",
            11,
            145,
        );
        svg.text(32.0, top + 40.0, &format!("Common numeric axis: {} per {unit}. Equal units do not establish matched exposure, loading, contrast or calibration.", report.left.signal_units), 11, 148);
        svg.text(32.0, top + 58.0, if report.left.source_sha256 != report.right.source_sha256 {
            "Different retained source identities; exposure matching is not established. Use CDF / QQ for shape comparison."
        } else {
            "Retained source identity is shared; loading/calibration comparability still needs review. Use CDF / QQ for shape."
        }, 11, 145);
    } else {
        let policy = if compatible {
            "Shared y scale"
        } else {
            "Independent L/R y scales; raw heights are NOT comparable"
        };
        svg.text(
            32.0,
            top + 24.0,
            &format!(
                "{policy} | L: {} | R: {} | NOT a density per {unit}",
                report.left.signal_units, report.right.signal_units
            ),
            12,
            135,
        );
    }
    let plot = Plot {
        x: 96.0,
        y: top + if density { 82.0 } else { 46.0 },
        w: 908.0,
        h: if density { 144.0 } else { 180.0 },
        lo: report.request.range[0],
        hi: report.request.range[1],
    };
    shade_gaps(svg, bins, plot);
    plot.axes(svg, ranges[0], &format!("Size ({unit})"));
    for i in 0..=4 {
        let f = f64::from(i) / 4.0;
        svg.label(
            plot.x + plot.w + 10.0,
            plot.y + (1.0 - f) * plot.h + 4.0,
            &number(ranges[1][0] + f * (ranges[1][1] - ranges[1][0])),
            11,
            10,
            RIGHT,
        );
    }
    svg.path("raw-left", &raw_path(bins, 0, plot, ranges[0]), LEFT, false);
    svg.path(
        "raw-right",
        &raw_path(bins, 1, plot, ranges[1]),
        RIGHT,
        true,
    );
    if density
        && bins
            .iter()
            .flat_map(|b| b.weights)
            .flatten()
            .all(|w| w == 0.0)
    {
        svg.text(32.0, top + 278.0, "All evaluated values are zero: axis 0..1 is a display fallback, not measured nonzero signal.", 10, 145);
    }
    Ok(())
}

/// Opt-in full comparison page. Legacy omission/explicit integrals retain the
/// exact old SVG. Invalid density geometry is refused, never clipped.
pub fn render_comparison_svg_with(
    report: &GelComparisonReport,
    mode: GelComparisonDisplay,
) -> Result<String, &'static str> {
    if mode == GelComparisonDisplay::SignalDensityPerSize {
        density_bins(report, &common_bins(report)?)?;
    }
    Ok(render_comparison_page(report, mode))
}

/// Profile panel only: GUI retains its native summaries and QQ instead of
/// duplicating the report page. Uses the same panel as the PDF projection.
pub fn render_comparison_profile_svg(
    report: &GelComparisonReport,
    mode: GelComparisonDisplay,
) -> Result<String, &'static str> {
    if report.left.size_kind != report.right.size_kind {
        return Err("Profile unavailable: incompatible size units.");
    }
    let bins = common_bins(report)?;
    let mut svg = Svg::new(340, "GENtle saved comparison profile");
    svg.text(32.0, 22.0, heading(mode), 17, 105);
    draw_profile_panel(&mut svg, report, &bins, mode, 22.0)?;
    svg.text(32.0, 314.0, "Teal solid = left | orange dashed = right | hatching = missing/excluded; no interpolation over gaps.", 11, 145);
    Ok(svg.finish())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gel_workbench::tests::{assert_xml, comparison, path};

    const DENSITY: GelComparisonDisplay = GelComparisonDisplay::SignalDensityPerSize;

    #[test]
    fn integral_mode_keeps_legacy_svg_bytes_and_report() {
        let report = comparison();
        let before = serde_json::to_value(&report).unwrap();
        assert_eq!(
            render_comparison_svg(&report),
            render_comparison_svg_with(&report, GelComparisonDisplay::BinIntegral).unwrap()
        );
        render_comparison_svg_with(&report, DENSITY).unwrap();
        assert_eq!(before, serde_json::to_value(report).unwrap());
    }

    #[test]
    fn unequal_split_and_offset_bins_preserve_flat_density_and_area() {
        // Invented integrated signals: every known interval has density 1.
        for edges in [
            vec![100.0, 125.0, 200.0, 400.0],
            vec![100.0, 200.0, 300.0, 400.0],
            vec![100.0, 125.0, 200.0, 300.0, 400.0],
        ] {
            let mut report = comparison();
            report.request.range = [100.0, 400.0];
            report.common_edges = edges;
            report.left_weights = report
                .common_edges
                .windows(2)
                .map(|w| Some(w[1] - w[0]))
                .collect();
            report.right_weights = report.left_weights.clone();
            report.excluded_intervals.clear();
            let bins = common_bins(&report).unwrap();
            let projected = density_bins(&report, &bins).unwrap();
            assert!(projected.iter().all(|b| b.weights == [Some(1.0); 2]));
            let area: f64 = projected
                .iter()
                .map(|b| b.weights[0].unwrap() * (b.end - b.start))
                .sum();
            assert!((area - 300.0).abs() < 1e-10);
            let svg = render_comparison_svg_with(&report, DENSITY).unwrap();
            assert_eq!(path(&svg, "raw-left"), path(&svg, "raw-right"));
            assert_eq!(
                path(&svg, "raw-left").matches(" 474.000").count(),
                2 * projected.len()
            );
            assert_xml(&svg);
        }
    }

    #[test]
    fn density_gaps_zero_axis_and_source_disclosures_are_explicit() {
        let mut report = comparison();
        report.right.source_sha256 = "b".repeat(64);
        let svg = render_comparison_profile_svg(&report, DENSITY).unwrap();
        assert!(svg.contains("per bp"));
        assert!(svg.contains("Common numeric axis"));
        assert!(svg.contains("Equal units do not establish matched exposure"));
        assert!(svg.contains("Different retained source identities"));
        assert!(svg.contains("exposure matching is not established"));
        assert!(!svg.contains("not exposure-matched"));
        assert_eq!(path(&svg, "raw-left").matches('M').count(), 3);
        assert!(path(&svg, "raw-left").ends_with("M822.400 248.000L1004.000 248.000"));
        assert!(svg.contains("fill=\"url(#missing)\""));
        assert!(!svg.contains("recorded distribution summary"));
        for profile in [&mut report.left, &mut report.right] {
            profile.size_kind = gentle_protocol::gel_image::GelSizeKind::SdsProteinKda;
        }
        assert!(
            render_comparison_profile_svg(&report, DENSITY)
                .unwrap()
                .contains("per kDa")
        );
        report.left_weights.fill(Some(0.0));
        report.right_weights.fill(Some(0.0));
        let svg = render_comparison_profile_svg(&report, DENSITY).unwrap();
        assert!(svg.contains("All evaluated values are zero"));
        assert!(svg.contains("display fallback, not measured nonzero signal"));
    }

    #[test]
    fn density_uses_common_numeric_scale_not_independent_maxima() {
        let mut report = comparison();
        report.right_weights = report
            .left_weights
            .iter()
            .map(|w| w.map(|w| w * 2.0))
            .collect();
        let svg = render_comparison_profile_svg(&report, DENSITY).unwrap();
        assert!(path(&svg, "raw-left").starts_with("M96.000 176.000"));
        assert!(path(&svg, "raw-right").starts_with("M96.000 104.000"));
        assert!(!svg.contains("Independent L/R"));
        report.right_weights.fill(Some(0.0));
        // A zero-only side must not force a 0..1 scale for small real signals.
        let svg = render_comparison_profile_svg(&report, DENSITY).unwrap();
        assert!(path(&svg, "raw-left").starts_with("M96.000 104.000"));
    }

    #[test]
    fn density_refuses_bad_values_or_geometry_without_clipping() {
        for value in [f64::INFINITY, f64::NAN, -1.0, f64::from_bits(1)] {
            let mut report = comparison();
            report.left_weights[0] = Some(value);
            assert!(render_comparison_profile_svg(&report, DENSITY).is_err());
            assert!(render_comparison_svg_with(&report, DENSITY).is_err());
        }
        let mut report = comparison();
        report.common_edges[1] = 100.0 + f64::EPSILON * 100.0;
        report.left_weights[0] = Some(f64::MAX);
        assert!(render_comparison_svg_with(&report, DENSITY).is_err());
        report = comparison();
        report.request.range = [-f64::MAX, f64::MAX];
        assert!(render_comparison_svg_with(&report, DENSITY).is_err());
    }
}
