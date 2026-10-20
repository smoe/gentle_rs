//! Pure, bounded SVG projections of frozen measured-gel reports.
//!
//! No source lookup, calibration, rebinning, quantile estimation or statistical
//! testing happens here. Cumulative geometry accumulates the supplied common-bin
//! integrals using the report's recorded totals. Missing support stays missing.

use gentle_protocol::gel_image::{GelBandSizingStatus, GelImageAnalysisReport};
use gentle_protocol::gel_workbench::{
    GelComparisonDisplay, GelComparisonReport, GelDistributionSummary, GelEvidenceKind,
    GelNumericSourceSummary, GelProfile, GelSourceSummary,
};
use std::fmt::Write;

const PAGE_WIDTH: usize = 1100;
const COMPARISON_HEIGHT: usize = 1460;
const MAX_BINS: usize = 20_000;
const MAX_QQ_POINTS: usize = 512;
const MAX_OVERVIEW_ITEMS: usize = 8;
const MAX_PREVIEW_BASE64: usize = 24 * 1024 * 1024;
const OVERVIEW_TOP: usize = 118;
const OVERVIEW_ROW: usize = 388;
const MAX_SOURCE_ITEMS: usize = 8;
const SOURCES_TOP: usize = 110;
const SOURCES_ROW: usize = 190;
const INK: &str = "#23313d";
const LEFT: &str = "#076d7b";
const RIGHT: &str = "#ad4c13";

mod algebra;
mod comparison_display;
pub use algebra::{render_algebra_provenance_svgs, render_algebra_svg};
pub use comparison_display::{render_comparison_profile_svg, render_comparison_svg_with};

fn xml(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            '\t' | '\n' | '\r' => out.push(ch),
            ch if ch >= ' ' && ch != '\u{fffe}' && ch != '\u{ffff}' => out.push(ch),
            _ => out.push('\u{fffd}'),
        }
    }
    out
}

fn short(value: &str, limit: usize) -> String {
    let mut chars = value.chars();
    let text: String = chars.by_ref().take(limit).collect();
    let mut text: String = text.split_whitespace().collect::<Vec<_>>().join(" ");
    if chars.next().is_some() {
        text.truncate(
            text.char_indices()
                .nth(limit.saturating_sub(3))
                .map_or(text.len(), |(i, _)| i),
        );
        text.push_str("...");
    }
    text
}

fn number(value: f64) -> String {
    if !value.is_finite() {
        "unavailable".into()
    } else if value == 0.0 {
        "0".into()
    } else if value.abs() < 0.01 || value.abs() >= 100_000.0 {
        format!("{value:.3e}")
    } else {
        format!("{value:.3}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .into()
    }
}

struct Svg(String, usize);

impl Svg {
    fn new(height: usize, title: &str) -> Self {
        let mut svg = Self(
            format!(
                "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{PAGE_WIDTH}\" height=\"{height}\" viewBox=\"0 0 {PAGE_WIDTH} {height}\" font-family=\"sans-serif\" fill=\"{INK}\" role=\"img\"><title>{}</title>\n\
             <defs><pattern id=\"missing\" width=\"8\" height=\"8\" patternUnits=\"userSpaceOnUse\"><rect width=\"8\" height=\"8\" fill=\"#f3f4f4\"/><path d=\"M0 8L8 0\" stroke=\"#d9dfe1\" stroke-width=\"1\"/></pattern></defs>\n",
                xml(title)
            ),
            0,
        );
        svg.rect(0.0, 0.0, PAGE_WIDTH as f64, height as f64, "#ffffff");
        svg
    }

    fn text(&mut self, x: f64, y: f64, value: &str, size: usize, limit: usize) {
        self.label(x, y, value, size, limit, INK);
    }

    fn label(&mut self, x: f64, y: f64, value: &str, size: usize, limit: usize, color: &str) {
        // Character limits alone cannot bound wide glyphs or unusual fonts.
        // Clip only visible text; the title retains the full supplied value.
        let width = (limit as f64 * size as f64 * 0.6)
            .min(PAGE_WIDTH as f64 - 32.0 - x)
            .max(0.0);
        let clip = self.1;
        self.1 += 1;
        let _ = writeln!(
            self.0,
            "<defs><clipPath id=\"text-clip-{clip}\"><rect x=\"{x:.3}\" y=\"{:.3}\" width=\"{width:.3}\" height=\"{:.3}\"/></clipPath></defs>",
            y - size as f64,
            size as f64 * 1.5
        );
        let _ = writeln!(
            self.0,
            "<text x=\"{x:.3}\" y=\"{y:.3}\" font-size=\"{size}\" fill=\"{color}\" clip-path=\"url(#text-clip-{clip})\"><title>{}</title>{}</text>",
            xml(value),
            xml(&short(value, limit))
        );
    }

    fn rect(&mut self, x: f64, y: f64, w: f64, h: f64, fill: &str) {
        let _ = writeln!(
            self.0,
            "<rect x=\"{x:.3}\" y=\"{y:.3}\" width=\"{w:.3}\" height=\"{h:.3}\" fill=\"{fill}\"/>"
        );
    }

    fn line(&mut self, x1: f64, y1: f64, x2: f64, y2: f64, color: &str) {
        let _ = writeln!(
            self.0,
            "<line x1=\"{x1:.3}\" y1=\"{y1:.3}\" x2=\"{x2:.3}\" y2=\"{y2:.3}\" stroke=\"{color}\"/>"
        );
    }

    fn path(&mut self, id: &str, data: &str, color: &str, dashed: bool) {
        let dash = if dashed {
            " stroke-dasharray=\"6 3\""
        } else {
            ""
        };
        let _ = writeln!(
            self.0,
            "<path id=\"{id}\" d=\"{data}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"2.3\"{dash}/>"
        );
    }

    fn finish(mut self) -> String {
        self.0.push_str("</svg>\n");
        self.0
    }
}

#[derive(Clone, Copy)]
struct Plot {
    x: f64,
    y: f64,
    w: f64,
    h: f64,
    lo: f64,
    hi: f64,
}

impl Plot {
    fn x(self, size: f64) -> f64 {
        self.x + (size - self.lo) / (self.hi - self.lo) * self.w
    }

    fn y(self, value: f64, range: [f64; 2]) -> f64 {
        self.y + (1.0 - (value - range[0]) / (range[1] - range[0])) * self.h
    }

    fn axes(self, svg: &mut Svg, y_range: [f64; 2], x_label: &str) {
        for i in 0..=4 {
            let fraction = f64::from(i) / 4.0;
            let y = self.y + (1.0 - fraction) * self.h;
            svg.line(self.x, y, self.x + self.w, y, "#e2e7e8");
            svg.label(
                self.x - 58.0,
                y + 4.0,
                &number(y_range[0] + fraction * (y_range[1] - y_range[0])),
                11,
                9,
                LEFT,
            );
            let x = self.x + fraction * self.w;
            svg.line(x, self.y + self.h, x, self.y + self.h + 5.0, INK);
            svg.text(
                x - 18.0,
                self.y + self.h + 22.0,
                &number(self.lo + fraction * (self.hi - self.lo)),
                11,
                11,
            );
        }
        svg.line(self.x, self.y, self.x, self.y + self.h, INK);
        svg.line(
            self.x,
            self.y + self.h,
            self.x + self.w,
            self.y + self.h,
            INK,
        );
        svg.text(
            self.x + self.w / 2.0 - 30.0,
            self.y + self.h + 44.0,
            x_label,
            12,
            24,
        );
    }
}

struct Bin {
    start: f64,
    end: f64,
    weights: [Option<f64>; 2],
}

fn common_bins(report: &GelComparisonReport) -> Result<Vec<Bin>, &'static str> {
    let n = report.left_weights.len();
    if n == 0 || n > MAX_BINS || report.excluded_intervals.len() > MAX_BINS {
        return Err(
            "Profile plots unavailable: empty data or render budget exceeded (20000 bins/intervals).",
        );
    }
    if report.right_weights.len() != n || report.common_edges.len() != n + 1 {
        return Err("Profile plots unavailable: common edges and weight counts do not match.");
    }
    let [lo, hi] = report.request.range;
    if !lo.is_finite()
        || !hi.is_finite()
        || hi <= lo
        || !(hi - lo).is_finite()
        || report
            .common_edges
            .iter()
            .any(|edge| !edge.is_finite() || *edge < lo || *edge > hi)
        || report
            .common_edges
            .windows(2)
            .any(|pair| pair[0] >= pair[1])
        || report
            .left_weights
            .iter()
            .chain(&report.right_weights)
            .flatten()
            .any(|v| !v.is_finite())
        || report
            .excluded_intervals
            .iter()
            .any(|[a, b]| !a.is_finite() || !b.is_finite() || a >= b)
    {
        return Err("Profile plots unavailable: invalid size range, bin geometry or weights.");
    }
    let mut excluded = report.excluded_intervals.clone();
    excluded.sort_by(|a, b| a[0].total_cmp(&b[0]));
    let mut gap = 0;
    Ok(report
        .common_edges
        .windows(2)
        .enumerate()
        .map(|(i, pair)| {
            while gap < excluded.len() && excluded[gap][1] <= pair[0] {
                gap += 1;
            }
            let masked = gap < excluded.len() && excluded[gap][0] < pair[1];
            Bin {
                start: pair[0],
                end: pair[1],
                weights: if masked {
                    [None, None]
                } else {
                    [report.left_weights[i], report.right_weights[i]]
                },
            }
        })
        .collect())
}

fn evidence(kind: GelEvidenceKind) -> &'static str {
    match kind {
        GelEvidenceKind::ImageDisplay => "Image-display evidence (not a vendor measurement)",
        GelEvidenceKind::VendorTrace => "Vendor-trace measurement (supplied measure/units)",
        GelEvidenceKind::SignedResidual => "Signed residual (derived; not a probability measure)",
        GelEvidenceKind::Synthetic => "Synthetic integrated bins (not experimental evidence)",
    }
}

fn profile_card(svg: &mut Svg, profile: &GelProfile, x: f64, side: &str, color: &str) {
    svg.rect(x, 101.0, 510.0, 181.0, "#f5f7f7");
    svg.rect(x, 101.0, 4.0, 181.0, color);
    svg.label(
        x + 12.0,
        125.0,
        &format!("{side}: {} | profile {}", profile.label, profile.id),
        15,
        54,
        color,
    );
    svg.text(
        x + 12.0,
        147.0,
        &format!("Run: {} | Lane: {}", profile.run_id, profile.lane_id),
        12,
        65,
    );
    svg.text(x + 12.0, 169.0, evidence(profile.evidence_kind), 12, 65);
    svg.text(
        x + 12.0,
        191.0,
        &format!(
            "Signal units: {} | size: {}",
            profile.signal_units,
            profile.size_kind.unit()
        ),
        12,
        65,
    );
    svg.text(
        x + 12.0,
        213.0,
        &format!("Measure: {}", profile.measure),
        12,
        65,
    );
    svg.text(
        x + 12.0,
        235.0,
        &format!("Source SHA-256: {}", profile.source_sha256),
        11,
        70,
    );
    let processing = profile
        .processing
        .first()
        .map_or("not specified", String::as_str);
    svg.text(
        x + 12.0,
        257.0,
        &format!(
            "Processing ({} steps): {processing}",
            profile.processing.len()
        ),
        11,
        70,
    );
}

fn weight_range(bins: &[Bin], side: usize) -> Option<[f64; 2]> {
    let mut lo: f64 = 0.0;
    let mut hi: f64 = 0.0;
    for value in bins.iter().filter_map(|bin| bin.weights[side]) {
        lo = lo.min(value);
        hi = hi.max(value);
    }
    if lo == hi {
        hi = 1.0;
    }
    (hi - lo).is_finite().then_some([lo, hi])
}

fn raw_path(bins: &[Bin], side: usize, plot: Plot, range: [f64; 2]) -> String {
    let mut path = String::new();
    // Separate horizontal subpaths preserve each integral, even for unequal
    // bin widths. No interpolation or line bridges over missing bins.
    for bin in bins {
        if let Some(weight) = bin.weights[side] {
            let y = plot.y(weight, range);
            let _ = write!(
                path,
                "M{:.3} {y:.3}L{:.3} {y:.3}",
                plot.x(bin.start),
                plot.x(bin.end)
            );
        }
    }
    path
}

fn cumulative_path(bins: &[Bin], side: usize, plot: Plot, total: f64) -> Option<String> {
    if !total.is_finite() || total <= 0.0 {
        return None;
    }
    let mut fraction = 0.0;
    let mut path = String::new();
    for bin in bins {
        if let Some(weight) = bin.weights[side] {
            if weight < 0.0 {
                return None;
            }
            let next = fraction + weight / total;
            if !next.is_finite() || next > 1.0 + 1e-8 {
                return None;
            }
            let _ = write!(
                path,
                "M{:.3} {:.3}L{:.3} {:.3}",
                plot.x(bin.start),
                plot.y(fraction.min(1.0), [0.0, 1.0]),
                plot.x(bin.end),
                plot.y(next.min(1.0), [0.0, 1.0])
            );
            fraction = next;
        }
    }
    // Refuse an inconsistent frozen normalizer rather than silently fitting it.
    ((fraction - 1.0).abs() <= 1e-8).then_some(path)
}

fn shade_gaps(svg: &mut Svg, bins: &[Bin], plot: Plot) {
    for bin in bins {
        if bin.weights.iter().any(Option::is_none) {
            svg.rect(
                plot.x(bin.start),
                plot.y,
                plot.x(bin.end) - plot.x(bin.start),
                plot.h,
                "url(#missing)",
            );
        }
    }
    if let Some(first) = bins.first().filter(|b| b.start > plot.lo) {
        svg.rect(
            plot.x,
            plot.y,
            plot.x(first.start) - plot.x,
            plot.h,
            "url(#missing)",
        );
    }
    if let Some(last) = bins.last().filter(|b| b.end < plot.hi) {
        svg.rect(
            plot.x(last.end),
            plot.y,
            plot.x + plot.w - plot.x(last.end),
            plot.h,
            "url(#missing)",
        );
    }
}

fn summary_card(
    svg: &mut Svg,
    profile: &GelProfile,
    summary: &GelDistributionSummary,
    threshold: Option<f64>,
    x: f64,
    color: &str,
) {
    svg.rect(x, 1100.0, 510.0, 112.0, "#f5f7f7");
    svg.label(
        x + 12.0,
        1122.0,
        &format!("{} | recorded distribution summary", profile.label),
        13,
        60,
        color,
    );
    svg.text(
        x + 12.0,
        1144.0,
        &format!(
            "Total weight: {} | {}",
            number(summary.total_weight),
            profile.signal_units
        ),
        12,
        65,
    );
    svg.text(
        x + 12.0,
        1166.0,
        &format!(
            "Q10: {} | median: {} | Q90: {} {}",
            number(summary.q10),
            number(summary.median),
            number(summary.q90),
            profile.size_kind.unit()
        ),
        12,
        65,
    );
    let fraction = match (threshold, summary.fraction_below_threshold) {
        (Some(t), Some(f)) => format!(
            "Weight fraction below {} {}: {}",
            number(t),
            profile.size_kind.unit(),
            number(f)
        ),
        (Some(t), None) => format!(
            "Weight fraction below {} {}: unavailable",
            number(t),
            profile.size_kind.unit()
        ),
        (None, _) => "Short-fragment threshold: not requested".into(),
    };
    svg.text(x + 12.0, 1188.0, &fraction, 12, 65);
}

fn warning_rows<'a>(
    svg: &mut Svg,
    rows: impl Iterator<Item = (&'a str, &'a str)>,
    x: f64,
    y: f64,
    count: usize,
    limit: usize,
) {
    let mut rows = rows;
    for index in 0..count {
        let Some((source, warning)) = rows.next() else {
            break;
        };
        svg.text(
            x,
            y + index as f64 * 18.0,
            &format!("{source}: {warning}"),
            11,
            limit,
        );
    }
    if rows.next().is_some() {
        svg.text(
            x,
            y + count as f64 * 18.0,
            "Additional warnings retained in the report; not shown on this page.",
            11,
            limit,
        );
    }
}

/// Render one fixed 1100 x 1460 page, using frozen common edges/weights and QQ
/// points only. Raw y values are per-bin integrated weights, never per-bp
/// densities. Incompatible raw measures get explicitly independent y scales.
/// Plot budgets fail visibly instead of decimating or inventing measurements.
/// Visible labels are shortened with full XML-safe values in hover titles.
pub fn render_comparison_svg(report: &GelComparisonReport) -> String {
    render_comparison_page(report, GelComparisonDisplay::BinIntegral)
}

fn render_comparison_page(report: &GelComparisonReport, mode: GelComparisonDisplay) -> String {
    let mut svg = Svg::new(COMPARISON_HEIGHT, "GENtle measured-gel comparison");
    let unit = report.left.size_kind.unit();
    svg.text(32.0, 34.0, "GENtle | Measured gel comparison", 23, 75);
    svg.text(
        32.0,
        60.0,
        &format!("Comparison: {}", report.request.comparison_id),
        14,
        105,
    );
    svg.text(
        32.0,
        83.0,
        &format!(
            "Requested interval: {} .. {} {unit} | frozen report: {}",
            number(report.request.range[0]),
            number(report.request.range[1]),
            report.schema
        ),
        12,
        130,
    );
    profile_card(&mut svg, &report.left, 32.0, "Left / solid", LEFT);
    profile_card(&mut svg, &report.right, 558.0, "Right / dashed", RIGHT);
    let distance_unit = if report.left.size_kind == report.right.size_kind {
        unit
    } else {
        "(incompatible size units)"
    };
    svg.text(
        32.0,
        312.0,
        &format!(
            "Descriptive D = {} (dimensionless) | Wasserstein distance = {} {distance_unit}",
            number(report.ks_distance),
            number(report.wasserstein_distance)
        ),
        17,
        105,
    );
    svg.text(32.0, 337.0, "p-value: not computed. Bins/pixels are correlated, not independent sampled fragments or biological replicates.", 13, 140);
    svg.text(32.0, 358.0, "Shape descriptors only: not a claim of yield, concentration, molecular identity or statistical significance.", 12, 140);
    svg.text(32.0, 392.0, comparison_display::heading(mode), 17, 105);
    svg.text(32.0, 700.0, "Normalized cumulative weight", 17, 60);
    svg.text(
        32.0,
        722.0,
        "Frozen totals on shared valid coverage; gaps stay gaps.",
        12,
        68,
    );
    svg.text(
        620.0,
        700.0,
        "QQ | stored quantiles, p = 0.01..0.99",
        16,
        52,
    );
    svg.text(
        620.0,
        722.0,
        "Same size units and scale on both axes; diagonal = equality.",
        11,
        68,
    );
    let mut notices: Vec<(&str, &str)> = Vec::new();
    let bins = common_bins(report);
    if report.left.size_kind != report.right.size_kind {
        notices.push(("Renderer", "Incompatible size units: all plots omitted; recorded metrics are not a valid same-unit comparison."));
        svg.text(
            96.0,
            520.0,
            "Plots unavailable: left and right sizes use incompatible units.",
            13,
            100,
        );
    } else if let Ok(bins) = &bins {
        let [lo, hi] = report.request.range;
        let cdf = Plot {
            x: 96.0,
            y: 750.0,
            w: 430.0,
            h: 280.0,
            lo,
            hi,
        };
        if let Err(message) =
            comparison_display::draw_profile_panel(&mut svg, report, bins, mode, 392.0)
        {
            notices.push(("Renderer", message));
        }
        shade_gaps(&mut svg, bins, cdf);
        cdf.axes(&mut svg, [0.0, 1.0], &format!("Size ({unit})"));
        for (side, total, color, kind, id) in [
            (
                0,
                report.left_summary.total_weight,
                LEFT,
                report.left.evidence_kind,
                "cdf-left",
            ),
            (
                1,
                report.right_summary.total_weight,
                RIGHT,
                report.right.evidence_kind,
                "cdf-right",
            ),
        ] {
            if let Some(path) = cumulative_path(bins, side, cdf, total)
                .filter(|_| kind != GelEvidenceKind::SignedResidual)
            {
                svg.path(id, &path, color, side == 1);
            } else {
                notices.push((
                    if side == 0 {
                        "Left cumulative"
                    } else {
                        "Right cumulative"
                    },
                    "Unavailable: signed measure, nonpositive total or inconsistent recorded normalization.",
                ));
            }
        }
        render_qq(
            &mut svg,
            report,
            Plot {
                x: 740.0,
                y: 750.0,
                w: 280.0,
                h: 280.0,
                lo,
                hi,
            },
        );
        svg.text(
            32.0,
            1090.0,
            "Linear within each known bin; no interpolation across gaps.",
            11,
            74,
        );
    } else if let Err(message) = bins {
        notices.push(("Renderer", message));
        svg.text(96.0, 520.0, message, 13, 115);
        svg.text(
            96.0,
            850.0,
            "Cumulative / QQ plots unavailable; inspect the retained report.",
            12,
            90,
        );
    }
    svg.text(32.0, 684.0, "Teal solid = left | orange dashed = right | hatching = missing/excluded in either profile; each curve keeps its own gaps.", 11, 145);
    summary_card(
        &mut svg,
        &report.left,
        &report.left_summary,
        report.request.short_threshold,
        32.0,
        LEFT,
    );
    summary_card(
        &mut svg,
        &report.right,
        &report.right_summary,
        report.request.short_threshold,
        558.0,
        RIGHT,
    );
    svg.text(32.0, 1240.0, "Evidence limits and warnings", 15, 90);
    if report.p_value.is_some() {
        notices.push(("Renderer", "Stored p-value ignored: correlated-bin comparisons do not support an inferential p-value."));
    }
    warning_rows(
        &mut svg,
        notices
            .into_iter()
            .chain(report.warnings.iter().map(|w| ("Comparison", w.as_str())))
            .chain(report.left.warnings.iter().map(|w| ("Left", w.as_str())))
            .chain(report.right.warnings.iter().map(|w| ("Right", w.as_str()))),
        32.0,
        1262.0,
        6,
        145,
    );
    svg.text(32.0, 1430.0, "Frozen inputs only. No re-extraction, recalibration, rebinning or inferential test at export. Rounded display values.", 11, 145);
    svg.finish()
}

fn render_qq(svg: &mut Svg, report: &GelComparisonReport, plot: Plot) {
    plot.axes(
        svg,
        [plot.lo, plot.hi],
        &format!("Left quantile ({})", report.left.size_kind.unit()),
    );
    svg.path(
        "qq-diagonal",
        &format!(
            "M{:.3} {:.3}L{:.3} {:.3}",
            plot.x,
            plot.y + plot.h,
            plot.x + plot.w,
            plot.y
        ),
        "#7c898f",
        true,
    );
    let _ = writeln!(
        svg.0,
        "<text transform=\"translate({:.3} {:.3}) rotate(-90)\" font-size=\"12\">Right quantile ({})</text>",
        plot.x - 66.0,
        plot.y + plot.h - 40.0,
        report.right.size_kind.unit()
    );
    if report.qq.len() > MAX_QQ_POINTS {
        svg.text(
            plot.x,
            plot.y + 90.0,
            "QQ omitted: more than 512 stored points.",
            12,
            38,
        );
        return;
    }
    let mut shown = 0;
    for point in &report.qq {
        if !(0.01..=0.99).contains(&point.probability)
            || !point.left.is_finite()
            || !point.right.is_finite()
            || !(plot.lo..=plot.hi).contains(&point.left)
            || !(plot.lo..=plot.hi).contains(&point.right)
        {
            continue;
        }
        let _ = writeln!(
            svg.0,
            "<circle class=\"qq-point\" cx=\"{:.3}\" cy=\"{:.3}\" r=\"2.5\" fill=\"{INK}\"><title>p={}; left={} {}; right={} {}</title></circle>",
            plot.x(point.left),
            plot.y(point.right, [plot.lo, plot.hi]),
            number(point.probability),
            number(point.left),
            report.left.size_kind.unit(),
            number(point.right),
            report.right.size_kind.unit()
        );
        shown += 1;
    }
    svg.text(
        620.0,
        1090.0,
        &format!(
            "{shown}/{} stored points shown; no quantiles recomputed.",
            report.qq.len()
        ),
        11,
        65,
    );
    if shown == 0 {
        svg.text(
            plot.x + 12.0,
            plot.y + 110.0,
            "No stored QQ points in this window.",
            12,
            36,
        );
    }
}

/// Render the explicitly ordered before/after items without inferring pairs
/// from lane names. The first two are side by side; up to four paired rows fit
/// on one bounded page (1100 x at most 1746). For more items the caller should
/// supply successive chunks of eight as separate overview pages; overflow is
/// disclosed, never squeezed into unreadable rows. Previews are embedded as
/// supplied, display-only, with their aspect ratio and orientation preserved.
/// Missing/invalid previews or base64 payloads over 24 MiB have visible placeholders.
/// Detailed analysis pages carry estimates/warnings beyond this short overview.
pub fn render_gel_overview_svg(title: &str, items: &[(&GelImageAnalysisReport, &str)]) -> String {
    let shown = items.len().min(MAX_OVERVIEW_ITEMS);
    let rows = shown.div_ceil(2);
    let height = OVERVIEW_TOP + rows.max(1) * OVERVIEW_ROW + 76;
    let mut svg = Svg::new(height, title);
    svg.text(32.0, 34.0, title, 23, 78);
    svg.text(
        32.0,
        60.0,
        "GENtle | Ordered before / after overview",
        15,
        100,
    );
    svg.text(32.0, 82.0, "Pairs follow supplied item order only, not lane-name matches or adjacency. No new measurements at export.", 12, 140);
    svg.text(32.0, 102.0, "Image-display evidence is not a vendor-trace measure. Band sizes are recorded calibration estimates only.", 12, 140);
    if shown == 0 {
        svg.text(
            32.0,
            160.0,
            "No analysis items supplied; no before/after pair inferred.",
            15,
            100,
        );
    }
    for (index, (report, preview)) in items.iter().take(MAX_OVERVIEW_ITEMS).enumerate() {
        let x = if index % 2 == 0 { 32.0 } else { 558.0 };
        let y = (OVERVIEW_TOP + index / 2 * OVERVIEW_ROW) as f64;
        overview_card(&mut svg, report, preview, index, x, y);
    }
    let footer = height as f64 - 48.0;
    let status = if items.len() > shown {
        format!(
            "Showing first {shown} of {} ordered items; remaining items require another overview page.",
            items.len()
        )
    } else if shown % 2 == 1 {
        format!(
            "{shown} ordered items. Final before slot has no supplied after partner; no missing image or result invented."
        )
    } else {
        format!(
            "{shown} ordered items | {} explicit before/after pairs. See detailed analysis and comparison pages.",
            shown / 2
        )
    };
    svg.text(32.0, footer, &status, 12, 140);
    svg.text(32.0, footer + 22.0, "Preview contrast does not establish abundance, recovery or molecular identity. Localization bounds are not confidence intervals.", 11, 145);
    svg.finish()
}

/// Render recorded original-document identities, not original document bytes.
/// At most eight ordered summaries fit on one 1100 x at most 1706 page. Both
/// 64-character SHA-256 values are visible in full-width monospace text, not
/// hover-only metadata. Other labels retain their full values in hover titles;
/// two warnings per source and a visible overflow notice keep rows bounded.
/// Callers can supply further chunks of eight as separate provenance pages.
/// This projection performs no source lookup, extraction or identity inference.
pub fn render_sources_svg(sources: &[GelSourceSummary]) -> String {
    let shown = sources.len().min(MAX_SOURCE_ITEMS);
    let height = SOURCES_TOP + shown.max(1) * SOURCES_ROW + 76;
    let mut svg = Svg::new(height, "GENtle | Original gel source provenance");
    svg.text(
        32.0,
        34.0,
        "GENtle | Original gel source provenance",
        23,
        78,
    );
    svg.text(32.0, 60.0, "Recorded PDF/native document identities and their linked extracted image IDs. No metadata inferred.", 13, 130);
    svg.text(32.0, 82.0, "Human review required: hashes bind content, not authenticated origin or scientific approval.", 12, 140);
    if shown == 0 {
        svg.text(
            32.0,
            160.0,
            "No original source summaries supplied; no document identity inferred.",
            15,
            105,
        );
    }
    for (index, source) in sources.iter().take(MAX_SOURCE_ITEMS).enumerate() {
        let y = (SOURCES_TOP + index * SOURCES_ROW) as f64;
        svg.rect(32.0, y, 1036.0, 176.0, "#f5f7f7");
        svg.rect(32.0, y, 4.0, 176.0, LEFT);
        svg.text(
            44.0,
            y + 22.0,
            &format!("{}. Source basename: {}", index + 1, source.source_name),
            15,
            108,
        );
        svg.text(
            44.0,
            y + 43.0,
            &format!("Linked image ID: {}", source.image_id),
            12,
            135,
        );
        svg.0.push_str("<g font-family=\"monospace\">\n");
        svg.text(
            44.0,
            y + 66.0,
            &format!("Original source SHA-256: {}", source.source_sha256),
            13,
            125,
        );
        svg.text(
            44.0,
            y + 87.0,
            &format!("Extracted image SHA-256: {}", source.extracted_image_sha256),
            13,
            125,
        );
        svg.0.push_str("</g>\n");
        svg.text(
            44.0,
            y + 108.0,
            &format!("Extraction: {}", source.extraction),
            12,
            135,
        );
        if source.warnings.is_empty() {
            svg.text(
                44.0,
                y + 130.0,
                "Warnings: none recorded; human review is still required.",
                11,
                145,
            );
        } else {
            warning_rows(
                &mut svg,
                source
                    .warnings
                    .iter()
                    .map(|warning| ("Source warning", warning.as_str())),
                44.0,
                y + 130.0,
                2,
                145,
            );
        }
    }
    let status = if sources.len() > shown {
        format!(
            "Showing first {shown} of {} sources; remaining sources require another provenance page.",
            sources.len()
        )
    } else {
        format!(
            "{shown} recorded source identities, in supplied order. Review extraction choices and source-specific warnings."
        )
    };
    svg.text(32.0, height as f64 - 48.0, &status, 12, 140);
    svg.text(
        32.0,
        height as f64 - 26.0,
        "Metadata only. Original document/archive payloads are not embedded on this page.",
        11,
        145,
    );
    svg.finish()
}

/// Render up to eight original CSV identities and their supplied profile
/// metadata on one bounded provenance page (1100 x at most 1706). SHA-256
/// values are printable in full-width monospace text; other labels are bounded
/// with full XML-safe hover titles. Callers paginate with chunks of eight.
/// No original payload, CSV parsing, source lookup or metadata inference occurs.
pub fn render_numeric_sources_svg(sources: &[GelNumericSourceSummary]) -> String {
    let shown = sources.len().min(MAX_SOURCE_ITEMS);
    let height = SOURCES_TOP + shown.max(1) * SOURCES_ROW + 76;
    let mut svg = Svg::new(height, "GENtle | Numeric gel source provenance");
    svg.text(32.0, 34.0, "GENtle | Numeric gel source provenance", 23, 78);
    svg.text(32.0, 60.0, "Recorded original CSV identities linked to frozen profile metadata. No signal data embedded.", 13, 130);
    svg.text(32.0, 82.0, "Human review required: hashes bind content, not authenticated origin or scientific approval.", 12, 140);
    if shown == 0 {
        svg.text(
            32.0,
            160.0,
            "No numeric source summaries supplied; no CSV identity inferred.",
            15,
            105,
        );
    }
    for (index, source) in sources.iter().take(MAX_SOURCE_ITEMS).enumerate() {
        let y = (SOURCES_TOP + index * SOURCES_ROW) as f64;
        svg.rect(32.0, y, 1036.0, 176.0, "#f5f7f7");
        svg.rect(32.0, y, 4.0, 176.0, LEFT);
        svg.text(
            44.0,
            y + 22.0,
            &format!("{}. Source basename: {}", index + 1, source.source_name),
            15,
            108,
        );
        svg.text(
            44.0,
            y + 43.0,
            &format!("Profile ID: {}", source.profile_id),
            12,
            135,
        );
        svg.0.push_str("<g font-family=\"monospace\">\n");
        svg.text(
            44.0,
            y + 66.0,
            &format!("Original CSV SHA-256: {}", source.sha256),
            13,
            125,
        );
        svg.0.push_str("</g>\n");
        svg.text(
            44.0,
            y + 87.0,
            &format!("Signal units: {}", source.signal_units),
            12,
            135,
        );
        svg.text(
            44.0,
            y + 108.0,
            &format!("Integration measure: {}", source.measure),
            12,
            135,
        );
    }
    let status = if sources.len() > shown {
        format!(
            "Showing first {shown} of {} numeric sources; remaining sources require another provenance page.",
            sources.len()
        )
    } else {
        format!(
            "{shown} recorded numeric source identities, in supplied order. Review declared units and integration measure."
        )
    };
    svg.text(32.0, height as f64 - 48.0, &status, 12, 140);
    svg.text(32.0, height as f64 - 26.0, "Metadata only. Original CSV payloads are not embedded on this page; no measurements are recomputed.", 11, 145);
    svg.finish()
}

fn safe_preview(preview: &str) -> bool {
    !preview.is_empty()
        && preview.len() <= MAX_PREVIEW_BASE64
        && preview.len() % 4 == 0
        && preview
            .trim_end_matches('=')
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'+' || b == b'/')
        && preview.len() - preview.trim_end_matches('=').len() <= 2
}

fn overview_card(
    svg: &mut Svg,
    report: &GelImageAnalysisReport,
    preview: &str,
    index: usize,
    x: f64,
    y: f64,
) {
    let color = if index % 2 == 0 { LEFT } else { RIGHT };
    let role = if index % 2 == 0 { "Before" } else { "After" };
    svg.rect(x, y, 510.0, 374.0, "#f5f7f7");
    svg.rect(x, y, 4.0, 374.0, color);
    svg.label(
        x + 12.0,
        y + 22.0,
        &format!(
            "{}. {role} | pair {} | {}",
            index + 1,
            index / 2 + 1,
            report.request.report_id
        ),
        14,
        57,
        color,
    );
    svg.text(
        x + 12.0,
        y + 43.0,
        &format!("Source: {}", report.image.source_name),
        12,
        65,
    );
    svg.text(
        x + 12.0,
        y + 61.0,
        &format!("SHA-256: {}", report.image.sha256),
        10,
        78,
    );
    svg.rect(x + 12.0, y + 73.0, 486.0, 126.0, "#e6eaea");
    if safe_preview(preview) {
        let _ = writeln!(
            svg.0,
            "<image x=\"{:.3}\" y=\"{:.3}\" width=\"486\" height=\"126\" preserveAspectRatio=\"xMidYMid meet\" href=\"data:image/png;base64,{}\"><title>Display-only preview of {}</title></image>",
            x + 12.0,
            y + 73.0,
            xml(preview),
            xml(&report.image.source_name)
        );
    } else {
        svg.text(
            x + 24.0,
            y + 138.0,
            "Preview unavailable: empty, invalid or over 24 MiB base64 budget.",
            12,
            62,
        );
    }
    svg.text(
        x + 12.0,
        y + 216.0,
        "Image-display evidence only; no vendor signal or abundance inferred.",
        11,
        74,
    );
    let calibrated = report
        .estimates
        .iter()
        .filter(|band| {
            band.status == GelBandSizingStatus::Interpolated
                && band
                    .estimated_size
                    .is_some_and(|size| size.is_finite() && size > 0.0)
        })
        .count();
    svg.text(
        x + 12.0,
        y + 238.0,
        &format!(
            "Calibrated: {calibrated}/{} bands ({}) | ladder: {}",
            report.estimates.len(),
            report.request.size_kind.unit(),
            report.request.ladder.label
        ),
        12,
        65,
    );
    for (i, band) in report.estimates.iter().take(3).enumerate() {
        let size = match (band.status, band.estimated_size) {
            (GelBandSizingStatus::Interpolated, Some(size)) if size.is_finite() && size > 0.0 => {
                format!("~{} {}", number(size), report.request.size_kind.unit())
            }
            (GelBandSizingStatus::OutsideCalibratedRange, _) => {
                "outside calibrated range; no extrapolation".into()
            }
            _ => "size unavailable".into(),
        };
        let bounds = band
            .localization_size_bounds
            .map(|[a, b]| {
                format!(
                    " | localization {}..{} {}",
                    number(a),
                    number(b),
                    report.request.size_kind.unit()
                )
            })
            .unwrap_or_default();
        svg.text(
            x + 12.0,
            y + 260.0 + i as f64 * 18.0,
            &format!("{size} | Lane {} | {}{bounds}", band.lane_id, band.label),
            11,
            74,
        );
    }
    if report.estimates.is_empty() {
        svg.text(
            x + 12.0,
            y + 260.0,
            "No reported sample-band estimates.",
            11,
            74,
        );
    } else if report.estimates.len() > 3 {
        svg.text(
            x + 12.0,
            y + 314.0,
            &format!(
                "+{} bands: see detailed analysis pages.",
                report.estimates.len() - 3
            ),
            11,
            74,
        );
    }
    warning_rows(
        svg,
        report
            .warnings
            .iter()
            .map(|w| ("Analysis warning", w.as_str()))
            .chain(
                report
                    .image
                    .warnings
                    .iter()
                    .map(|w| ("Image warning", w.as_str())),
            )
            .chain(
                report
                    .estimates
                    .iter()
                    .flat_map(|band| band.warnings.iter().map(|w| ("Band warning", w.as_str()))),
            ),
        x + 12.0,
        y + 337.0,
        1,
        74,
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use gentle_protocol::gel_image::{
        GelBandSizeEstimate, GelImageAnalysisRequest, GelImageDescriptor, GelImageLadder,
        GelImagePoint, GelMigrationDirection, GelSizeKind,
    };
    use gentle_protocol::gel_workbench::{
        GelCompareRequest, GelNumericSource, GelProfileFileRequest, GelQqPoint,
        GelSourceImportRequest, GelSourceRecord,
    };
    use svg::node::element::tag::Type;
    use svg::parser::Event;

    // Hand-crafted, synthetic reports only. These constructors are the complete
    // deterministic provenance; no external scientific fixtures are implied.
    fn profile(label: &str) -> GelProfile {
        GelProfile {
            id: label.into(),
            label: label.into(),
            run_id: "synthetic-run".into(),
            lane_id: "A1".into(),
            source_sha256: "a".repeat(64),
            evidence_kind: GelEvidenceKind::ImageDisplay,
            size_kind: GelSizeKind::LinearDnaBp,
            signal_units: "display intensity x pixels".into(),
            edges: vec![100.0, 600.0],
            weights: vec![Some(999.0)],
            measure: "synthetic integrated display signal".into(),
            processing: vec![],
            warnings: vec![],
        }
    }

    pub(super) fn comparison() -> GelComparisonReport {
        let summary = GelDistributionSummary {
            total_weight: 20.0,
            q10: 110.0,
            median: 200.0,
            q90: 580.0,
            fraction_below_threshold: Some(0.5),
        };
        GelComparisonReport {
            schema: "synthetic".into(),
            request: GelCompareRequest {
                comparison_id: "synthetic-comparison".into(),
                left_profile: "before".into(),
                right_profile: "after".into(),
                range: [100.0, 600.0],
                short_threshold: Some(300.0),
            },
            left: profile("before"),
            right: profile("after"),
            common_edges: vec![100.0, 200.0, 400.0, 500.0, 600.0],
            left_weights: vec![Some(10.0), None, Some(10.0), Some(0.0)],
            right_weights: vec![Some(5.0), None, Some(10.0), Some(5.0)],
            excluded_intervals: vec![[200.0, 400.0]],
            left_summary: summary.clone(),
            right_summary: summary,
            qq: vec![GelQqPoint {
                probability: 0.5,
                left: 300.0,
                right: 500.0,
            }],
            ks_distance: 0.25,
            wasserstein_distance: 42.0,
            p_value: None,
            warnings: vec![],
        }
    }

    fn analysis() -> GelImageAnalysisReport {
        GelImageAnalysisReport {
            schema: "synthetic".into(),
            algorithm: "synthetic".into(),
            gentle_version: "synthetic".into(),
            image: GelImageDescriptor {
                schema: "synthetic".into(),
                image_id: "image".into(),
                source_name: "synthetic.png".into(),
                sha256: "b".repeat(64),
                byte_count: 0,
                width: 200,
                height: 100,
                format: "png".into(),
                color_type: "rgb8".into(),
                tiff_page: None,
                warnings: vec!["synthetic display warning".into()],
            },
            request: GelImageAnalysisRequest {
                report_id: "synthetic-analysis".into(),
                image_id: "image".into(),
                image_sha256: "b".repeat(64),
                size_kind: GelSizeKind::LinearDnaBp,
                migration: GelMigrationDirection::Down,
                lanes: vec![],
                ladder: GelImageLadder {
                    lane_id: "L1".into(),
                    label: "synthetic ladder".into(),
                    source: "hand-crafted test".into(),
                    gel_system: None,
                    prestained: false,
                    bands: vec![],
                },
                sample_bands: vec![],
            },
            calibration_bands: vec![],
            estimates: vec![GelBandSizeEstimate {
                band_id: "sample".into(),
                lane_id: "A1".into(),
                label: "synthetic band".into(),
                center: GelImagePoint { x: 10.0, y: 20.0 },
                status: GelBandSizingStatus::Interpolated,
                estimated_size: Some(250.0),
                reference_band_ids: vec![],
                localization_size_bounds: Some([240.0, 260.0]),
                warnings: vec![],
            }],
            warnings: vec![],
        }
    }

    pub(super) fn assert_xml(document: &str) {
        let mut stack = Vec::new();
        for event in svg::read(document).expect("read SVG") {
            match event {
                Event::Error(error) => panic!("invalid SVG/XML: {error}"),
                Event::Tag(name, Type::Start, _) => stack.push(name),
                Event::Tag(name, Type::End, _) => assert_eq!(stack.pop(), Some(name)),
                _ => {}
            }
        }
        assert!(stack.is_empty());
        assert!(!document.contains('\0'));
    }

    pub(super) fn path<'a>(document: &'a str, id: &str) -> &'a str {
        document
            .split_once(&format!("<path id=\"{id}\" d=\""))
            .expect("path exists")
            .1
            .split_once('"')
            .unwrap()
            .0
    }

    #[test]
    fn legacy_comparison_svg_bytes() {
        // FNV-1a is a portable byte-regression checksum, not an evidence digest.
        let checksum = render_comparison_svg(&comparison())
            .bytes()
            .fold(0xcbf29ce484222325_u64, |hash, byte| {
                (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3)
            });
        assert_eq!(checksum, 9_151_467_362_457_725_946);
    }

    #[test]
    fn raw_uses_common_integrals_not_original_bins_or_density() {
        let mut report = comparison();
        report.left_weights = vec![Some(10.0), Some(10.0), None, Some(0.0)];
        report.excluded_intervals.clear();
        let document = render_comparison_svg(&report);
        // The 100-wide and 200-wide bins have the same height, not 10/width.
        assert!(
            path(&document, "raw-left")
                .starts_with("M96.000 438.000L277.600 438.000M277.600 438.000L640.800 438.000")
        );
        assert!(!document.contains("999"));
        assert!(document.contains("integrated weight per common bin"));
        assert_xml(&document);
    }

    #[test]
    fn raw_and_cumulative_do_not_bridge_missing_or_excluded_bins() {
        let document = render_comparison_svg(&comparison());
        assert_eq!(
            path(&document, "raw-left"),
            "M96.000 438.000L277.600 438.000M640.800 438.000L822.400 438.000M822.400 618.000L1004.000 618.000"
        );
        assert_eq!(
            path(&document, "cdf-left"),
            "M96.000 1030.000L182.000 890.000M354.000 890.000L440.000 750.000M440.000 750.000L526.000 750.000"
        );
        assert!(document.contains("fill=\"url(#missing)\""));
        let mut report = comparison();
        report.left_weights[1] = Some(800.0);
        report.right_weights[1] = Some(900.0);
        assert_eq!(document, render_comparison_svg(&report));
    }

    #[test]
    fn qq_is_same_unit_square_and_only_supplied_probabilities() {
        let mut report = comparison();
        for probability in [0.0, 0.009, 0.01, 0.99, 0.991, 1.0] {
            report.qq.push(GelQqPoint {
                probability,
                left: 300.0,
                right: 500.0,
            });
        }
        let document = render_comparison_svg(&report);
        assert_eq!(document.matches("class=\"qq-point\"").count(), 3);
        assert!(document.contains("cx=\"852.000\" cy=\"806.000\""));
        assert_eq!(
            path(&document, "qq-diagonal"),
            "M740.000 1030.000L1020.000 750.000"
        );
        assert!(document.contains("Left quantile (bp)"));
        assert!(document.contains("Right quantile (bp)"));
    }

    #[test]
    fn evidence_units_and_non_inference_are_explicit() {
        let mut report = comparison();
        report.right.evidence_kind = GelEvidenceKind::VendorTrace;
        report.right.signal_units = "FU x seconds".into();
        report.right.measure = "vendor integrated time trace".into();
        report.p_value = Some(0.000001);
        let document = render_comparison_svg(&report);
        for phrase in [
            "Image-display evidence (not a vendor measurement)",
            "Vendor-trace measurement",
            "FU x seconds",
            "vendor integrated time trace",
            "Independent L/R y scales",
            "p-value: not computed",
            "Bins/pixels are correlated",
            "D = 0.25",
            "Wasserstein distance = 42 bp",
            "Stored p-value ignored",
            "Run: synthetic-run | Lane: A1",
        ] {
            assert!(document.contains(phrase), "missing {phrase}");
        }
        assert!(!document.contains("1.000e-6"));
    }

    #[test]
    fn user_strings_are_xml_safe_and_full_labels_survive_in_hover() {
        let malicious = "<&\"'\0";
        let mut report = comparison();
        report.request.comparison_id = malicious.into();
        report.left.label = "long biological label ".repeat(30);
        report.left.run_id = malicious.into();
        report.left.lane_id = malicious.into();
        report.left.signal_units = malicious.into();
        report.left.measure = malicious.into();
        report.left.source_sha256 = malicious.into();
        report.left.processing = vec![malicious.into()];
        report.warnings = vec![malicious.into()];
        let document = render_comparison_svg(&report);
        assert!(document.contains("&lt;&amp;&quot;&apos;\u{fffd}"));
        assert!(document.contains(&format!(
            "<title>Left / solid: {} | profile before</title>",
            report.left.label
        )));
        assert!(!document.contains(malicious));
        assert_xml(&document);
        let mut image = analysis();
        image.image.source_name = malicious.into();
        image.image.sha256 = malicious.into();
        image.estimates[0].label = malicious.into();
        image.warnings = vec![malicious.into()];
        let overview = render_gel_overview_svg(malicious, &[(&image, "AAAA")]);
        assert!(overview.contains("&lt;&amp;&quot;&apos;\u{fffd}"));
        assert_xml(&overview);
    }

    #[test]
    fn invalid_numeric_geometry_and_budgets_fail_visibly() {
        let mut report = comparison();
        report.common_edges[1] = f64::NAN;
        let document = render_comparison_svg(&report);
        assert!(document.contains("invalid size range, bin geometry or weights"));
        assert!(!document.contains("NaN"));
        assert!(!document.contains("id=\"raw-left\""));
        report = comparison();
        report.left_weights = vec![Some(1.0); MAX_BINS + 1];
        assert!(render_comparison_svg(&report).contains("render budget exceeded"));
        report = comparison();
        report.left_summary.total_weight = 99.0;
        let document = render_comparison_svg(&report);
        assert!(!document.contains("id=\"cdf-left\""));
        assert!(document.contains("inconsistent recorded normalization"));
        assert!(document.contains("Total weight: 99"));
    }

    #[test]
    fn full_engine_bin_budget_is_rendered_without_decimation() {
        let mut report = comparison();
        report.request.range = [100.0, 100.0 + MAX_BINS as f64];
        report.common_edges = (0..=MAX_BINS).map(|i| 100.0 + i as f64).collect();
        report.left_weights = vec![Some(1.0); MAX_BINS];
        report.right_weights = report.left_weights.clone();
        report.excluded_intervals.clear();
        report.left_summary.total_weight = MAX_BINS as f64;
        report.right_summary.total_weight = MAX_BINS as f64;
        let document = render_comparison_svg(&report);
        assert_eq!(MAX_BINS, 20_000);
        assert_eq!(path(&document, "raw-left").matches('M').count(), MAX_BINS);
        assert_eq!(path(&document, "raw-right").matches('M').count(), MAX_BINS);
        assert_eq!(path(&document, "cdf-left").matches('M').count(), MAX_BINS);
        assert!(!document.contains("render budget exceeded"));
    }

    #[test]
    fn preview_budget_accepts_large_verified_sources() {
        let mut preview = "A".repeat(MAX_PREVIEW_BASE64);
        assert_eq!(MAX_PREVIEW_BASE64, 24 * 1024 * 1024);
        assert!(safe_preview(&preview));
        preview.push_str("AAAA");
        assert!(!safe_preview(&preview));
        assert!(!safe_preview("AAAA\" onload=\"alert(1)"));
    }

    #[test]
    fn overview_preserves_order_original_preview_and_bounded_pairs() {
        let mut reports = vec![analysis(); 10];
        for (i, report) in reports.iter_mut().enumerate() {
            report.request.report_id = format!("ordered-report-{i}");
        }
        let items: Vec<_> = reports.iter().map(|report| (report, "AAAA")).collect();
        let document = render_gel_overview_svg("Synthetic overview", &items);
        assert!(document.contains("width=\"1100\" height=\"1746\""));
        assert_eq!(document.matches("<image ").count(), 8);
        assert_eq!(
            document
                .matches("href=\"data:image/png;base64,AAAA\"")
                .count(),
            8
        );
        assert!(document.contains("preserveAspectRatio=\"xMidYMid meet\""));
        assert!(document.contains("1. Before | pair 1 | ordered-report-0"));
        assert!(document.contains("2. After | pair 1 | ordered-report-1"));
        assert!(document.contains("3. Before | pair 2 | ordered-report-2"));
        assert!(
            document.find("ordered-report-0").unwrap() < document.find("ordered-report-1").unwrap()
        );
        assert!(!document.contains("ordered-report-8"));
        assert!(document.contains("Showing first 8 of 10"));
        assert!(document.contains("~250 bp"));
        assert!(document.contains("localization 240..260 bp"));
        assert!(document.contains("synthetic display warning"));
        assert!(document.contains(&"b".repeat(64)));
        assert_xml(&document);
    }

    #[test]
    fn overview_empty_unpaired_and_overflow_content_are_disclosed() {
        let empty = render_gel_overview_svg("Empty", &[]);
        assert!(empty.contains("No analysis items supplied"));
        assert_xml(&empty);
        let mut report = analysis();
        report.estimates = vec![report.estimates[0].clone(); 100];
        report.warnings = vec!["warning ".repeat(100); 100];
        let document = render_gel_overview_svg("Odd", &[(&report, "\" onload=\"attack")]);
        assert!(!document.contains("<image "));
        assert!(document.contains("Preview unavailable"));
        assert!(document.contains("no supplied after partner"));
        assert!(document.contains("+97 bands"));
        assert!(document.contains("Additional warnings retained"));
        assert!(document.contains("height=\"582\""));
        assert_xml(&document);
    }

    #[test]
    fn outside_calibration_is_not_a_fabricated_size_and_protein_keeps_units() {
        let mut report = analysis();
        report.request.size_kind = GelSizeKind::SdsProteinKda;
        report.estimates[0].status = GelBandSizingStatus::OutsideCalibratedRange;
        let document = render_gel_overview_svg("Protein", &[(&report, "AAAA")]);
        assert!(document.contains("Calibrated: 0/1 bands (kDa)"));
        assert!(document.contains("outside calibrated range; no extrapolation"));
        assert!(!document.contains("~250"));
    }

    #[test]
    fn source_provenance_prints_full_hashes_without_original_archive_payload() {
        let source = GelSourceRecord {
            image_id: "synthetic-extracted-image".into(),
            source_name: "synthetic-original.D5000".into(),
            source_sha256: "0123456789abcdef".repeat(4),
            original_base64: "UEsDBHN5bnRoZXRpYy1wcml2YXRlLWFyY2hpdmUtcGF5bG9hZA==".into(),
            extracted_image_sha256: "fedcba9876543210".repeat(4),
            extraction: "explicit display-only native thumbnail".into(),
            warnings: vec!["Synthetic warning: review thumbnail contrast.".into()],
            request: GelSourceImportRequest {
                image_id: "synthetic-extracted-image".into(),
                path: "synthetic-original.D5000".into(),
                pdf_object: None,
                allow_native_thumbnail: true,
            },
        };
        let document = render_sources_svg(&[GelSourceSummary::from(&source)]);
        assert!(document.contains("synthetic-original.D5000"));
        assert!(document.contains("Linked image ID: synthetic-extracted-image"));
        assert!(document.contains("explicit display-only native thumbnail"));
        assert!(document.contains("review thumbnail contrast"));
        assert!(document.contains("Human review required"));
        assert!(document.contains("font-family=\"monospace\""));
        for (label, hash) in [
            ("Original source SHA-256", &source.source_sha256),
            ("Extracted image SHA-256", &source.extracted_image_sha256),
        ] {
            assert_eq!(hash.len(), 64);
            // Assert visible text, not just the full value retained in a title.
            assert!(document.contains(&format!("</title>{label}: {hash}</text>")));
        }
        assert!(!document.contains(&source.original_base64));
        assert!(!document.contains("synthetic-private-archive-payload"));
        assert!(!document.contains("original_base64"));
        assert!(!document.contains("base64"));
        assert!(!document.contains("href="));
        assert!(!document.contains("<image"));
        assert_xml(&document);
    }

    #[test]
    fn source_provenance_escapes_metadata_and_bounds_rows_and_warnings() {
        let source = GelSourceSummary {
            image_id: "<&\"'\0".into(),
            source_name: "synthetic<&.pdf".into(),
            source_sha256: "a".repeat(64),
            extracted_image_sha256: "b".repeat(64),
            extraction: "PDF raster <object> & explicit selection".into(),
            warnings: vec!["Review <source> & extraction ".repeat(100); 20],
        };
        let mut sources = vec![source; 9];
        sources[8].source_name = "must-not-appear.pdf".into();
        let document = render_sources_svg(&sources);
        assert!(document.contains("width=\"1100\" height=\"1706\""));
        assert!(document.contains("&lt;&amp;&quot;&apos;\u{fffd}"));
        assert!(document.contains("synthetic&lt;&amp;.pdf"));
        assert!(document.contains("PDF raster &lt;object&gt; &amp; explicit selection"));
        assert!(document.contains("Showing first 8 of 9 sources"));
        assert!(!document.contains("must-not-appear.pdf"));
        assert_eq!(
            document
                .matches("</title>Additional warnings retained")
                .count(),
            8
        );
        assert!(document.contains(&format!(
            "<title>Source warning: {}</title>",
            xml(&sources[0].warnings[0])
        )));
        assert_xml(&document);
        let empty = render_sources_svg(&[]);
        assert!(empty.contains("No original source summaries supplied"));
        assert_xml(&empty);
        sources.truncate(1);
        sources[0].warnings.clear();
        assert!(
            render_sources_svg(&sources)
                .contains("Warnings: none recorded; human review is still required.")
        );
    }

    #[test]
    fn numeric_source_provenance_prints_full_hash_without_csv_payload() {
        let source = GelNumericSource {
            source_name: "synthetic-original.csv".into(),
            sha256: "0123456789abcdef".repeat(4),
            original_base64: "c3ludGhldGljLXByaXZhdGUtY3N2LXBheWxvYWQ=".into(),
            request: GelProfileFileRequest {
                path: "synthetic-original.csv".into(),
                profile_id: "synthetic-numeric-profile".into(),
                label: "synthetic trace".into(),
                run_id: "synthetic-run".into(),
                lane_id: "A1".into(),
                size_kind: GelSizeKind::LinearDnaBp,
                signal_units: "FU x seconds".into(),
                measure: "synthetic integrated time trace".into(),
            },
        };
        let summary = GelNumericSourceSummary {
            profile_id: source.request.profile_id.clone(),
            source_name: source.source_name.clone(),
            sha256: source.sha256.clone(),
            signal_units: source.request.signal_units.clone(),
            measure: source.request.measure.clone(),
        };
        let document = render_numeric_sources_svg(&[summary]);
        assert!(document.contains("synthetic-original.csv"));
        assert!(document.contains("Profile ID: synthetic-numeric-profile"));
        assert!(document.contains("Signal units: FU x seconds"));
        assert!(document.contains("Integration measure: synthetic integrated time trace"));
        assert!(document.contains("font-family=\"monospace\""));
        assert_eq!(source.sha256.len(), 64);
        assert!(document.contains(&format!(
            "</title>Original CSV SHA-256: {}</text>",
            source.sha256
        )));
        assert!(!document.contains(&source.original_base64));
        assert!(!document.contains("synthetic-private-csv-payload"));
        assert!(!document.contains("original_base64"));
        assert!(!document.contains("base64"));
        assert!(!document.contains("href="));
        assert!(!document.contains("<image"));
        assert_xml(&document);
    }

    #[test]
    fn numeric_source_provenance_escapes_metadata_and_bounds_ordered_pages() {
        let source = GelNumericSourceSummary {
            profile_id: "<&\"'\0".into(),
            source_name: "synthetic<&.csv".into(),
            sha256: "a".repeat(64),
            signal_units: "FU & ticks".into(),
            measure: "Explicit <integral> & declared measure ".repeat(100),
        };
        let mut sources = vec![source; 9];
        sources[7].profile_id = "eighth-profile".into();
        sources[8].source_name = "must-not-appear.csv".into();
        let document = render_numeric_sources_svg(&sources);
        assert!(document.contains("width=\"1100\" height=\"1706\""));
        assert!(document.contains("&lt;&amp;&quot;&apos;\u{fffd}"));
        assert!(document.contains("synthetic&lt;&amp;.csv"));
        assert!(document.contains("Signal units: FU &amp; ticks"));
        assert!(document.contains(&format!(
            "<title>Integration measure: {}</title>",
            xml(&sources[0].measure)
        )));
        assert!(document.contains("Showing first 8 of 9 numeric sources"));
        assert!(!document.contains("must-not-appear.csv"));
        assert!(
            document.find("synthetic&lt;&amp;.csv").unwrap()
                < document.find("eighth-profile").unwrap()
        );
        assert_eq!(
            document.matches("</title>Original CSV SHA-256: ").count(),
            8
        );
        assert_xml(&document);
        let empty = render_numeric_sources_svg(&[]);
        assert!(empty.contains("No numeric source summaries supplied"));
        assert_xml(&empty);
    }
}
