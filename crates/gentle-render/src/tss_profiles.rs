//! Readable, deterministic SVG pages from accession-pinned TSS profile reports.
//!
//! No sequence lookup, scoring, calibration or correlation is performed here.
//! Row scales use the report's raw values, with clipping applied only to display.
//! All horizontal geometry uses transcript-oriented motif-window starts. PFM
//! logos use uncorrected Shannon information against a uniform DNA background,
//! not the scoring model's pseudocounts and never a consensus approximation.
//! Source-provided selection labels and non-claims repeat with selected panels;
//! neither selection evidence nor missing annotation releases are inferred.

use gentle_protocol::tss_profiles::{
    NON_CLAIMS, REPORT_SCHEMA, ResolvedTssMatrix, TssCalibrationState, TssGeometry,
    TssMatrixComparison, TssPeak, TssProfileRenderOptions, TssProfileReport, TssProfileTrack,
    TssProfileWindow, TssScaleMode, TssStrand,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;
use std::ops::Range;
use svg::Node;
use svg::node::element::{Circle, Group, Line, Path, Rectangle, Text, Title};

// Match the canonical locus-evidence page and its shared genomic plot frame.
// This keeps detailed TSS pages horizontally registered with the context page
// when both are viewed as one paginated report.
const PAGE_WIDTH: f64 = 1400.0;
const MAX_PAGE_HEIGHT: f64 = 11_200.0;
const MAX_PANELS_PER_PAGE: usize = 32;
const MARGIN: f64 = 34.0;
const TEXT_WIDTH: f64 = PAGE_WIDTH - 2.0 * MARGIN;
const LABEL_WIDTH: f64 = 205.0;
const PLOT_LEFT: f64 = 255.0;
const PLOT_RIGHT: f64 = 1050.0;
const PLOT_WIDTH: f64 = PLOT_RIGHT - PLOT_LEFT;
const PLOT_TOP: f64 = 24.0;
const PLOT_HEIGHT: f64 = 128.0;
const MIN_ROW_HEIGHT: f64 = 224.0;
const ROW_GAP: f64 = 16.0;
const PANEL_GAP: f64 = 28.0;
const LOGO_COLUMNS: usize = 20;
const LOGO_COLUMN_WIDTH: f64 = 8.5;
const LOGO_HEIGHT: f64 = 64.0;
const LOGO_SEGMENT_HEIGHT: f64 = 90.0;
const COLORS: [&str; 6] = [
    "#176b87", "#a34e24", "#30794d", "#8b4160", "#846416", "#485f91",
];

/// One self-contained SVG page. Page numbers/counts are one-based and per gene.
/// Promoter IDs follow their first appearance on this page; a continued TSS can
/// appear on multiple pages. Rendering never reorders the report's windows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TssRenderedPage {
    pub gene_id: String,
    pub gene_symbol: String,
    pub page_number: usize,
    pub page_count: usize,
    pub promoter_ids: Vec<String>,
    pub svg: String,
}

/// Project a verified report into independent SVG documents without mutating it.
///
/// The effective scale is the option override or the panel's typed policy
/// (normally independent per matrix and TSS). Cross-matrix shared ranges require
/// `CrossSourceCalibrated`, one exact nonempty ID, lowercase 64-hex SHA-256 and a
/// nonempty statement. This checks the supplied binding, not its scientific merit.
///
/// `panels_per_page` must be 1..=32 and is an upper bound, not a request to shrink
/// rows. Pages are 1400 pixels wide and at most 11,200 pixels tall (62.7 MB for a
/// single RGBA raster at native size; higher raster scales require their own
/// budget). Oversized TSSs continue at whole-row/comparison boundaries, repeating
/// their identity and axes. A single indivisible block that cannot fit is an error.
/// Invalid geometry, score arrays, PFMs or ambiguous matrix identities also fail
/// before any pages are returned. Empty verified window lists yield no pages.
pub fn render_tss_profile_pages(
    report: &TssProfileReport,
    options: &TssProfileRenderOptions,
) -> Result<Vec<TssRenderedPage>, String> {
    let scale = validate_report(report, options)?;
    let footer = footer_block(report, scale);
    let mut genes: Vec<Vec<&TssProfileWindow>> = Vec::new();
    let mut gene_indices = BTreeMap::new();
    for window in &report.windows {
        let index = *gene_indices
            .entry(window.record.gene_id.as_str())
            .or_insert_with(|| {
                genes.push(Vec::new());
                genes.len() - 1
            });
        genes[index].push(window);
    }
    let mut output = Vec::new();
    for windows in genes {
        let record = &windows[0].record;
        let title = TextBlock::new(
            &format!("{} | {}", record.gene_symbol, record.gene_id),
            TEXT_WIDTH,
            22.0,
        );
        let header = page_header(report, 1, 1);
        let page_top = MARGIN + title.height() + 12.0 + header.height() + PANEL_GAP;
        let fixed_height = page_top + PANEL_GAP + footer.height() + MARGIN;
        let capacity = MAX_PAGE_HEIGHT - fixed_height;
        if capacity <= 0.0 {
            return Err("TSS page header/footer exceeds the 10000-pixel height budget".into());
        }
        let layouts = windows
            .iter()
            .enumerate()
            .map(|(index, window)| {
                WindowLayout::new(report, window, scale, index + 1, windows.len())
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut fragments = Vec::new();
        for (index, layout) in layouts.iter().enumerate() {
            fragments.extend(layout.fragments(index, capacity)?);
        }
        let mut pages: Vec<Vec<Fragment>> = Vec::new();
        let mut used = 0.0;
        for fragment in fragments {
            let height = layouts[fragment.window].fragment_height(&fragment);
            let needs_page = pages.last().is_none_or(|page| {
                page.len() >= options.panels_per_page || used + PANEL_GAP + height > capacity
            });
            if needs_page {
                pages.push(Vec::new());
                used = 0.0;
            }
            let page = pages.last_mut().expect("a page was just allocated");
            if !page.is_empty() {
                used += PANEL_GAP;
            }
            used += height;
            page.push(fragment);
        }
        let page_count = pages.len();
        for (page_index, fragments) in pages.into_iter().enumerate() {
            let content_height = fragments
                .iter()
                .map(|fragment| layouts[fragment.window].fragment_height(fragment))
                .sum::<f64>()
                + PANEL_GAP * fragments.len().saturating_sub(1) as f64;
            let height = (fixed_height + content_height).ceil();
            let mut body = Group::new();
            body.append(
                Rectangle::new()
                    .set("width", PAGE_WIDTH)
                    .set("height", height)
                    .set("fill", "#ffffff"),
            );
            title.draw(&mut body, MARGIN, MARGIN, "gene-title");
            page_header(report, page_index + 1, page_count).draw(
                &mut body,
                MARGIN,
                MARGIN + title.height() + 12.0,
                "page-header",
            );
            let mut y = page_top;
            let mut promoter_ids = Vec::new();
            for fragment in &fragments {
                let layout = &layouts[fragment.window];
                let id = &layout.window.record.promoter_id;
                if !promoter_ids.contains(id) {
                    promoter_ids.push(id.clone());
                }
                body.append(layout.draw(fragment, y));
                y += layout.fragment_height(fragment) + PANEL_GAP;
            }
            footer.draw(&mut body, MARGIN, y, "calibration-nonclaims-footer");
            let svg = svg::Document::new()
                .set("width", PAGE_WIDTH)
                .set("height", height)
                .set("viewBox", (0, 0, PAGE_WIDTH, height))
                .set("font-family", "monospace")
                .set("font-size", 14)
                .set("fill", "#233441")
                .set("role", "img")
                .set("aria-labelledby", "tss-page-title")
                .set(
                    "data-score-kind",
                    report.panel_resolution.panel.score_kind.as_str(),
                )
                .set("data-scale-mode", scale_name(scale))
                .set("data-page-number", page_index + 1)
                .set("data-page-count", page_count)
                .set("data-gentle-plot-left", PLOT_LEFT)
                .set("data-gentle-plot-right", PLOT_RIGHT)
                .add(
                    Title::new(format!(
                        "{} TSS profiles, page {} of {}",
                        record.gene_symbol,
                        page_index + 1,
                        page_count
                    ))
                    .set("id", "tss-page-title"),
                )
                .add(body)
                .to_string();
            output.push(TssRenderedPage {
                gene_id: record.gene_id.clone(),
                gene_symbol: record.gene_symbol.clone(),
                page_number: page_index + 1,
                page_count,
                promoter_ids,
                svg,
            });
        }
    }
    Ok(output)
}

fn validate_report(
    report: &TssProfileReport,
    options: &TssProfileRenderOptions,
) -> Result<TssScaleMode, String> {
    if !(1..=MAX_PANELS_PER_PAGE).contains(&options.panels_per_page) {
        return Err("panels_per_page must be between 1 and 32".into());
    }
    if report.schema != REPORT_SCHEMA {
        return Err(format!("Unsupported TSS report schema: {}", report.schema));
    }
    let resolution = &report.panel_resolution;
    let panel = &resolution.panel;
    let scale = options.scale_mode.unwrap_or(panel.scale_mode);
    score_units(&panel.score_kind)?;
    let cross_matrix = scale == TssScaleMode::Shared && resolution.matrices.len() > 1;
    if cross_matrix && panel.calibration_state != TssCalibrationState::CrossSourceCalibrated {
        return Err("Shared cross-matrix scales require typed CrossSourceCalibrated evidence; matching units or prose are insufficient".into());
    }
    if panel.calibration_state == TssCalibrationState::CrossSourceCalibrated {
        let id_ok = panel
            .calibration_id
            .as_deref()
            .is_some_and(|id| !id.trim().is_empty() && id == id.trim());
        let digest_ok = panel.calibration_sha256.as_deref().is_some_and(|digest| {
            digest.len() == 64
                && digest
                    .bytes()
                    .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        });
        if !id_ok || !digest_ok || panel.calibration_statement.trim().is_empty() {
            return Err("CrossSourceCalibrated requires an exact nonempty calibration_id, lowercase 64-hex calibration_sha256 and nonempty statement".into());
        }
    }
    if resolution.matrices.is_empty() && !report.windows.is_empty() {
        return Err("TSS report contains windows but no resolved matrices".into());
    }
    let mut matrices = BTreeMap::new();
    let mut orders = BTreeSet::new();
    for matrix in &resolution.matrices {
        let spec = &matrix.specification;
        if spec.source_id.is_empty()
            || matrices.insert(spec.source_id.as_str(), matrix).is_some()
            || !orders.insert(spec.display_order)
        {
            return Err(format!(
                "Ambiguous resolved matrix identity/order: {}",
                spec.source_id
            ));
        }
        if spec
            .score_kind
            .as_deref()
            .is_some_and(|kind| kind != panel.score_kind)
        {
            return Err(format!(
                "Matrix {} has a conflicting score_kind",
                spec.source_id
            ));
        }
        if matrix.matrix_counts.is_empty()
            || matrix.matrix_counts.iter().any(|column| {
                column
                    .iter()
                    .any(|value| !value.is_finite() || *value < 0.0)
                    || !column.iter().any(|value| *value > 0.0)
            })
        {
            return Err(format!(
                "Matrix {} requires a finite, nonnegative full PFM with nonempty columns; consensus is not a logo",
                spec.source_id
            ));
        }
        matrix_color(matrix, 0)?;
    }
    let mut promoters = BTreeSet::new();
    let mut symbols = BTreeMap::new();
    for window in &report.windows {
        let record = &window.record;
        let context = &record.promoter_id;
        if record.gene_id.is_empty() || context.is_empty() || !promoters.insert(context) {
            return Err(format!(
                "Missing gene identity or duplicate promoter: {context}"
            ));
        }
        if symbols
            .insert(&record.gene_id, &record.gene_symbol)
            .is_some_and(|symbol| symbol != &record.gene_symbol)
        {
            return Err(format!("Conflicting gene symbols for {}", record.gene_id));
        }
        let geometry = &record.geometry;
        let length = geometry
            .length()
            .ok_or_else(|| format!("Window length overflow: {context}"))?;
        let genomic_length = geometry
            .end_1based
            .checked_sub(geometry.start_1based)
            .and_then(|n| n.checked_add(1));
        if geometry.start_1based == 0
            || u64::try_from(length).ok() != genomic_length
            || geometry.genomic_at(geometry.upstream_bp) != Some(geometry.tss_1based)
            || geometry.relative_at(length - 1).is_none()
        {
            return Err(format!(
                "Inconsistent transcript/genomic geometry: {context}"
            ));
        }
        let mut seen = BTreeSet::new();
        for track in &window.tracks {
            let matrix = matrices
                .get(track.accession.as_str())
                .ok_or_else(|| format!("Unresolved track {} in {context}", track.accession))?;
            let valid_length = length
                .checked_sub(track.motif_length_bp)
                .map_or(0, |n| n + 1);
            if !seen.insert(&track.accession)
                || track.motif_length_bp != matrix.matrix_counts.len()
                || track.forward_scores.len() != valid_length
                || track.reverse_scores.len() != valid_length
            {
                return Err(format!(
                    "Duplicate track or inconsistent motif/window-start lengths for {} in {context}",
                    track.accession
                ));
            }
            if track
                .forward_scores
                .iter()
                .chain(&track.reverse_scores)
                .flatten()
                .any(|score| !score.is_finite())
            {
                return Err(format!(
                    "Non-finite raw score for {} in {context}; unavailable scores must be null",
                    track.accession
                ));
            }
            for peak in track
                .forward_maximum
                .iter()
                .chain(track.reverse_maximum.iter())
                .chain(&track.forward_peaks)
                .chain(&track.reverse_peaks)
            {
                if !peak.score.is_finite() || peak.local_start_0based >= valid_length {
                    return Err(format!(
                        "Invalid reported peak for {} in {context}",
                        track.accession
                    ));
                }
            }
        }
        if seen.len() != matrices.len() {
            return Err(format!("Missing resolved matrix tracks in {context}"));
        }
        for comparison in &window.comparisons {
            for accession in [&comparison.left_accession, &comparison.right_accession] {
                if matrices
                    .get(accession.as_str())
                    .is_none_or(|matrix| matrix.specification.factor_id != comparison.factor_id)
                {
                    return Err(format!(
                        "Comparison {accession} does not belong to declared factor {}",
                        comparison.factor_id
                    ));
                }
            }
            if comparison
                .pearson
                .iter()
                .chain(comparison.spearman.iter())
                .any(|value| !value.is_finite())
            {
                return Err("Non-finite comparison; undefined metrics must be null".into());
            }
        }
    }
    Ok(scale)
}

fn score_units(kind: &str) -> Result<&'static str, String> {
    match kind {
        "llr_bits" => Ok("bits (matrix-specific LLR)"),
        "true_log_odds_bits" => Ok("bits (true log odds)"),
        "llr_quantile" | "true_log_odds_quantile" => Ok("unitless matrix-specific quantile"),
        "llr_background_quantile" | "true_log_odds_background_quantile" => {
            Ok("unitless modeled background quantile")
        }
        "llr_background_tail_log10" | "true_log_odds_background_tail_log10" => {
            Ok("-log10(modeled background-tail probability)")
        }
        _ => Err(format!("Unsupported TSS score_kind: {kind}")),
    }
}

fn scale_name(scale: TssScaleMode) -> &'static str {
    match scale {
        TssScaleMode::Independent => "independent",
        TssScaleMode::Shared => "shared",
    }
}

fn matrix_color(matrix: &ResolvedTssMatrix, index: usize) -> Result<String, String> {
    let Some(color) = matrix.specification.color_hint.as_deref() else {
        return Ok(COLORS[index % COLORS.len()].into());
    };
    let safe = color.strip_prefix('#').is_some_and(|hex| {
        matches!(hex.len(), 3 | 4 | 6 | 8) && hex.bytes().all(|byte| byte.is_ascii_hexdigit())
    }) || (!color.is_empty() && color.bytes().all(|byte| byte.is_ascii_alphabetic()));
    if !safe {
        return Err(format!(
            "Unsupported color_hint for {}: use a hex or named color, never a resource URL",
            matrix.specification.source_id
        ));
    }
    Ok(color.into())
}

fn page_header(report: &TssProfileReport, number: usize, count: usize) -> TextBlock {
    let resolution = &report.panel_resolution;
    TextBlock::new(
        &format!(
            "TSS TFBS profiles | page {number} of {count}\nPanel {}: {} | {} resolved matrices\nReference: {} | assembly: {} | annotation: {}\nPanel SHA-256: {}",
            resolution.panel.panel_id,
            resolution.panel.label,
            resolution.matrices.len(),
            report.reference.genome_id,
            report.reference.assembly,
            report
                .reference
                .annotation_release
                .as_deref()
                .unwrap_or("not separately declared"),
            resolution.panel_sha256,
        ),
        TEXT_WIDTH,
        14.0,
    )
}

fn footer_block(report: &TssProfileReport, scale: TssScaleMode) -> TextBlock {
    let panel = &report.panel_resolution.panel;
    let mut text = String::from("CALIBRATION AND INTERPRETATION\n");
    text.push_str(match scale {
        TssScaleMode::Independent => "Independent matrix-specific ranges, recomputed per TSS. Equal trace heights across rows or TSS pages do not mean equal scores.\n",
        TssScaleMode::Shared => "Shared range within each TSS, recomputed for each TSS; not a common range across TSS pages.\n",
    });
    let _ = writeln!(
        text,
        "Calibration state: {} | ID: {} | SHA-256: {}",
        match panel.calibration_state {
            TssCalibrationState::MatrixSpecific => "matrix_specific",
            TssCalibrationState::CrossSourceCalibrated =>
                "cross_source_calibrated (declared report binding, not independently verified)",
        },
        panel.calibration_id.as_deref().unwrap_or("none"),
        panel.calibration_sha256.as_deref().unwrap_or("none")
    );
    let _ = writeln!(
        text,
        "Statement: {}",
        if panel.calibration_statement.trim().is_empty() {
            "No cross-matrix calibration statement supplied."
        } else {
            &panel.calibration_statement
        }
    );
    text.push_str(if panel.clip_negative {
        "Display: max(raw score, 0). Negative values are clipped only in the figure; raw scores and stored comparisons are unchanged.\n"
    } else {
        "Display: raw scores, including negatives. Stored comparisons are unchanged.\n"
    });
    text.push_str("PFM logos: p(base) * [2 - Shannon entropy] bits; uniform DNA reference, no pseudocount or small-sample correction. Logo bits are not track scores.\n");
    text.push_str("Full per-matrix normalization references are embedded in SVG row titles and retained in the source report.\n");
    for (key, value) in &report.score_policy {
        let _ = writeln!(text, "Score policy {key}: {value}");
    }
    let _ = writeln!(text, "Verification: {}", report.verification);
    for warning in &report.warnings {
        let _ = writeln!(text, "Warning: {warning}");
    }
    text.push_str(NON_CLAIMS);
    if !report.non_claims.trim().is_empty() && report.non_claims.trim() != NON_CLAIMS {
        let _ = write!(text, "\nReport interpretation: {}", report.non_claims);
    }
    TextBlock::new(&text, TEXT_WIDTH, 13.0)
}

#[derive(Debug)]
struct TextBlock {
    lines: Vec<String>,
    width: f64,
    size: f64,
}

fn text_width(value: &str, size: f64) -> f64 {
    value
        .chars()
        .map(|c| if c.is_ascii() { size * 0.64 } else { size })
        .sum()
}

impl TextBlock {
    fn new(value: &str, width: f64, size: f64) -> Self {
        let mut lines = Vec::new();
        for paragraph in value.split('\n') {
            let mut line = String::new();
            for word in paragraph.split_whitespace() {
                if !line.is_empty()
                    && text_width(&line, size) + text_width(word, size) + size * 0.64 > width
                {
                    lines.push(std::mem::take(&mut line));
                }
                if !line.is_empty() {
                    line.push(' ');
                }
                for character in word.chars() {
                    let advance = if character.is_ascii() {
                        size * 0.64
                    } else {
                        size
                    };
                    if !line.is_empty() && text_width(&line, size) + advance > width {
                        lines.push(std::mem::take(&mut line));
                    }
                    line.push(character);
                }
            }
            lines.push(line);
        }
        Self { lines, width, size }
    }

    fn height(&self) -> f64 {
        self.lines.len() as f64 * (self.size + 5.0)
    }

    fn draw(&self, parent: &mut Group, x: f64, y: f64, role: &str) {
        let mut block = Group::new()
            .set("data-role", role)
            .set("data-x", x)
            .set("data-y", y)
            .set("data-width", self.width)
            .set("data-height", self.height());
        for (index, line) in self.lines.iter().enumerate() {
            if !line.is_empty() {
                block.append(text_node(
                    x,
                    y + self.size + index as f64 * (self.size + 5.0),
                    line,
                    self.size,
                ));
            }
        }
        parent.append(block);
    }
}

fn text_node(x: f64, y: f64, value: &str, size: f64) -> Text {
    // Explicit advances bound text independently of the installed fallback font.
    Text::new(value)
        .set("x", x)
        .set("y", y)
        .set("font-size", size)
        .set("textLength", text_width(value, size))
        .set("lengthAdjust", "spacingAndGlyphs")
}

#[derive(Clone, Copy, Debug)]
struct ScoreRange {
    min: f64,
    max: f64,
}

impl ScoreRange {
    fn from_tracks<'a>(tracks: impl Iterator<Item = &'a TssProfileTrack>, clip: bool) -> Self {
        let mut range = Self { min: 0.0, max: 0.0 };
        for score in tracks
            .flat_map(|track| track.forward_scores.iter().chain(&track.reverse_scores))
            .flatten()
        {
            let score = display_score(*score, clip);
            range.min = range.min.min(score);
            range.max = range.max.max(score);
        }
        if range.min == range.max {
            range.max = 1.0;
        }
        range
    }

    fn y(self, score: f64, top: f64) -> f64 {
        // Normalize before subtraction so opposite extreme finite scores cannot overflow.
        let magnitude = self.min.abs().max(self.max.abs());
        let fraction = ((score / magnitude - self.min / magnitude)
            / (self.max / magnitude - self.min / magnitude))
            .clamp(0.0, 1.0);
        top + (1.0 - fraction) * PLOT_HEIGHT
    }

    fn ticks(self) -> [f64; 3] {
        [self.min, self.min * 0.5 + self.max * 0.5, self.max]
    }
}

fn display_score(score: f64, clip: bool) -> f64 {
    if clip { score.max(0.0) } else { score }
}

fn number(value: f64) -> String {
    if value == 0.0 {
        return "0".into();
    }
    if value.abs() >= 1_000_000.0 || value.abs() < 0.001 {
        return format!("{value:.3e}");
    }
    format!("{value:.4}")
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string()
}

#[derive(Clone, Copy)]
struct LocalAxis {
    length: usize,
}

impl LocalAxis {
    fn x(self, index: f64) -> f64 {
        if self.length == 1 {
            PLOT_LEFT + PLOT_WIDTH * 0.5
        } else {
            PLOT_LEFT + index / (self.length - 1) as f64 * PLOT_WIDTH
        }
    }

    fn interval(self, start: usize, end: usize) -> (f64, f64) {
        if self.length == 1 {
            return (PLOT_LEFT, PLOT_WIDTH);
        }
        let left = self.x((start as f64 - 0.5).max(0.0));
        let right = self.x((end as f64 - 0.5).min((self.length - 1) as f64));
        (left, (right - left).max(0.0))
    }
}

struct AxisTick {
    index: usize,
    relative: String,
    genomic: String,
    tier: usize,
    center: f64,
    width: f64,
}

struct AxisLayout {
    axis: LocalAxis,
    ticks: Vec<AxisTick>,
    tiers: usize,
}

impl AxisLayout {
    fn new(geometry: &TssGeometry) -> Self {
        let length = geometry.length().expect("validated geometry");
        let axis = LocalAxis { length };
        let widest_genomic = geometry
            .start_1based
            .to_string()
            .len()
            .max(geometry.end_1based.to_string().len()) as f64
            * 13.0
            * 0.64;
        let desired_step = ((length - 1) as f64 * (widest_genomic + 48.0) / PLOT_WIDTH).max(1.0);
        let magnitude = 10_f64.powf(desired_step.log10().floor());
        let step = ([1.0, 2.0, 5.0, 10.0]
            .into_iter()
            .find(|n| *n * magnitude >= desired_step)
            .unwrap_or(10.0)
            * magnitude) as i64;
        let step = step.max(1);
        let low = -(geometry.upstream_bp as i64);
        let high = geometry.downstream_bp as i64;
        let mut indices = BTreeSet::from([0, geometry.upstream_bp, length - 1]);
        let mut relative = low / step * step;
        while relative <= high {
            indices.insert((relative - low) as usize);
            let Some(next) = relative.checked_add(step) else {
                break;
            };
            relative = next;
        }
        let mut tier_ends = Vec::<f64>::new();
        let mut ticks = Vec::new();
        for index in indices {
            let relative = geometry
                .relative_at(index)
                .expect("validated relative coordinate");
            let relative = if relative > 0 {
                format!("+{relative}")
            } else {
                relative.to_string()
            };
            let genomic = geometry
                .genomic_at(index)
                .expect("validated genomic coordinate")
                .to_string();
            let width = text_width(&relative, 13.0).max(text_width(&genomic, 13.0));
            let center = axis.x(index as f64).clamp(
                PLOT_LEFT + width * 0.5,
                PLOT_LEFT + PLOT_WIDTH - width * 0.5,
            );
            let tier = tier_ends
                .iter()
                .position(|end| *end + 16.0 <= center - width * 0.5)
                .unwrap_or(tier_ends.len());
            if tier == tier_ends.len() {
                tier_ends.push(0.0);
            }
            tier_ends[tier] = center + width * 0.5;
            ticks.push(AxisTick {
                index,
                relative,
                genomic,
                tier,
                center,
                width,
            });
        }
        Self {
            axis,
            ticks,
            tiers: tier_ends.len(),
        }
    }

    fn height(&self) -> f64 {
        16.0 + self.tiers as f64 * 46.0
    }

    fn draw(&self, group: &mut Group, y: f64) {
        group.append(
            Line::new()
                .set("x1", PLOT_LEFT)
                .set("x2", PLOT_LEFT + PLOT_WIDTH)
                .set("y1", y)
                .set("y2", y)
                .set("stroke", "#8296a0"),
        );
        for tier in 0..self.tiers {
            group.append(text_node(
                MARGIN,
                y + 21.0 + tier as f64 * 46.0,
                "TSS-relative motif start (bp)",
                13.0,
            ));
            group.append(text_node(
                MARGIN,
                y + 40.0 + tier as f64 * 46.0,
                "Genomic base (1-based)",
                13.0,
            ));
        }
        for tick in &self.ticks {
            let x = self.axis.x(tick.index as f64);
            group.append(
                Line::new()
                    .set("x1", x)
                    .set("x2", x)
                    .set("y1", y)
                    .set("y2", y + 6.0)
                    .set("stroke", "#526b79"),
            );
            for (role, label, offset) in [
                ("relative-tick", &tick.relative, 21.0),
                ("genomic-tick", &tick.genomic, 40.0),
            ] {
                group.append(
                    text_node(
                        tick.center,
                        y + offset + tick.tier as f64 * 46.0,
                        label,
                        13.0,
                    )
                    .set("text-anchor", "middle")
                    .set("data-role", role)
                    .set("data-local-start", tick.index)
                    .set("data-value", label.as_str())
                    .set("data-axis-x", x)
                    .set("data-label-width", tick.width)
                    .set(
                        "font-weight",
                        if tick.relative == "0" {
                            "bold"
                        } else {
                            "normal"
                        },
                    ),
                );
            }
        }
    }
}

struct RowLayout<'a> {
    matrix: &'a ResolvedTssMatrix,
    track: &'a TssProfileTrack,
    label: TextBlock,
    identity: TextBlock,
    notes: TextBlock,
    range: ScoreRange,
    color: String,
    logo_top: f64,
    height: f64,
}

impl<'a> RowLayout<'a> {
    fn new(
        report: &TssProfileReport,
        window: &TssProfileWindow,
        matrix: &'a ResolvedTssMatrix,
        track: &'a TssProfileTrack,
        range: ScoreRange,
        index: usize,
    ) -> Result<Self, String> {
        let spec = &matrix.specification;
        let label = TextBlock::new(&spec.label, LABEL_WIDTH, 15.0);
        let identity = TextBlock::new(
            &format!(
                "{} | {}\n{} bp | {}\n{}",
                spec.source_id,
                spec.factor_id,
                track.motif_length_bp,
                report.panel_resolution.panel.score_kind,
                score_units(&report.panel_resolution.panel.score_kind)?,
            ),
            LABEL_WIDTH,
            12.0,
        );
        let valid = |scores: &[Option<f64>]| scores.iter().flatten().count();
        let terminal = window.record.geometry.length().expect("validated geometry")
            - track.forward_scores.len();
        let mut notes = format!(
            "Valid starts: local + {}/{}; local - {}/{}; gray terminal positions: {terminal} unscored.\nRaw maxima: local + {}; local - {}.",
            valid(&track.forward_scores),
            track.forward_scores.len(),
            valid(&track.reverse_scores),
            track.reverse_scores.len(),
            peak_label(track.forward_maximum.as_ref(), &window.record.geometry),
            peak_label(track.reverse_maximum.as_ref(), &window.record.geometry),
        );
        if track
            .forward_scores
            .iter()
            .chain(&track.reverse_scores)
            .flatten()
            .next()
            .is_none()
        {
            notes.push_str("\nNo valid scored windows; the 0..1 axis is a display fallback, not evidence of zero signal.");
        } else if track
            .forward_scores
            .iter()
            .chain(&track.reverse_scores)
            .flatten()
            .all(|value| display_score(*value, report.panel_resolution.panel.clip_negative) == 0.0)
        {
            notes.push_str(
                "\nAll valid displayed values are zero; missing windows remain unavailable.",
            );
        }
        if track.normalization_reference.is_null() {
            notes.push_str("\nNormalization: not supplied; no additional calibration inferred.");
        } else {
            let reference = &track.normalization_reference;
            let background = reference
                .get("background_model")
                .and_then(serde_json::Value::as_str)
                .unwrap_or("not specified");
            let model = reference
                .get("chance_model")
                .and_then(serde_json::Value::as_str)
                .unwrap_or("not specified");
            let _ = write!(notes, "\nBackground: {background} | model: {model}");
        }
        let notes = TextBlock::new(&notes, PLOT_WIDTH, 13.0);
        let logo_top = 16.0 + label.height() + 8.0 + identity.height() + 26.0;
        let left_height = logo_top
            + matrix.matrix_counts.len().div_ceil(LOGO_COLUMNS) as f64 * LOGO_SEGMENT_HEIGHT
            + 14.0;
        let right_height = PLOT_TOP + PLOT_HEIGHT + 24.0 + notes.height() + 14.0;
        Ok(Self {
            matrix,
            track,
            label,
            identity,
            notes,
            range,
            color: matrix_color(matrix, index)?,
            logo_top,
            height: MIN_ROW_HEIGHT.max(left_height).max(right_height),
        })
    }

    fn draw(&self, y: f64, axis: LocalAxis, geometry: &TssGeometry, clip: bool) -> Group {
        let mut row = Group::new()
            .set("data-role", "matrix-row")
            .set(
                "data-accession",
                self.matrix.specification.source_id.as_str(),
            )
            .set(
                "data-display-order",
                self.matrix.specification.display_order,
            )
            .set("data-matrix-sha256", self.matrix.matrix_sha256.as_str())
            .set("data-matrix-version", self.matrix.version.as_str())
            .set("data-y-min", self.range.min)
            .set("data-y-max", self.range.max)
            .set("data-y", y)
            .set("data-height", self.height);
        row.append(Title::new(format!(
            "{}; matrix SHA-256 {}; normalization reference (stored, not recomputed): {}",
            self.matrix.specification.source_id,
            self.matrix.matrix_sha256,
            self.track.normalization_reference,
        )));
        row.append(
            Rectangle::new()
                .set("x", MARGIN - 12.0)
                .set("y", y)
                .set("width", TEXT_WIDTH + 24.0)
                .set("height", self.height)
                .set("rx", 4)
                .set("fill", "#f8fafb"),
        );
        self.label.draw(&mut row, MARGIN, y + 16.0, "matrix-label");
        self.identity.draw(
            &mut row,
            MARGIN,
            y + 24.0 + self.label.height(),
            "matrix-identity-units",
        );
        row.append(text_node(
            MARGIN,
            y + self.logo_top - 9.0,
            "PFM information logo (bits)",
            13.0,
        ));
        draw_logo(&mut row, self.matrix, MARGIN, y + self.logo_top);
        let top = y + PLOT_TOP;
        row.append(
            Rectangle::new()
                .set("x", PLOT_LEFT)
                .set("y", top)
                .set("width", PLOT_WIDTH)
                .set("height", PLOT_HEIGHT)
                .set("fill", "#ffffff")
                .set("stroke", "#b4c3ca"),
        );
        let first_unscored = self.track.forward_scores.len();
        if first_unscored < axis.length {
            let (x, width) = axis.interval(first_unscored, axis.length);
            row.append(
                Rectangle::new()
                    .set("data-role", "terminal-unscored")
                    .set("data-local-start", first_unscored)
                    .set("data-local-end-exclusive", axis.length)
                    .set("x", x)
                    .set("y", top)
                    .set("width", width)
                    .set("height", PLOT_HEIGHT)
                    .set("fill", "#dfe6ea")
                    .add(Title::new(format!(
                        "Unscored base positions {}..{}; no complete {}-bp motif window",
                        first_unscored,
                        axis.length - 1,
                        self.track.motif_length_bp
                    ))),
            );
        }
        for (strand, scores) in [
            (TssStrand::Plus, &self.track.forward_scores),
            (TssStrand::Minus, &self.track.reverse_scores),
        ] {
            let mut index = 0;
            while index < scores.len() {
                if scores[index].is_some() {
                    index += 1;
                    continue;
                }
                let start = index;
                while index < scores.len() && scores[index].is_none() {
                    index += 1;
                }
                let (x, width) = axis.interval(start, index);
                row.append(Rectangle::new().set("data-role", "unavailable-windows")
                    .set("data-local-strand", strand.as_str()).set("data-local-start", start)
                    .set("data-local-end-exclusive", index)
                    .set("x", x).set("y", top).set("width", width).set("height", PLOT_HEIGHT)
                    .set("fill", "#eacb80").set("fill-opacity", 0.26)
                    .add(Title::new(format!("Unavailable motif-local {} windows (including ambiguous bases), not zero scores", strand.as_str()))));
            }
        }
        for value in self.range.ticks() {
            let yy = self.range.y(value, top);
            row.append(
                Line::new()
                    .set("x1", PLOT_LEFT)
                    .set("x2", PLOT_LEFT + PLOT_WIDTH)
                    .set("y1", yy)
                    .set("y2", yy)
                    .set("stroke", "#d5dfe4")
                    .set("stroke-width", 0.8),
            );
            row.append(
                text_node(PLOT_RIGHT + 8.0, yy + 4.0, &number(value), 13.0)
                    .set("text-anchor", "start")
                    .set("data-role", "score-tick")
                    .set("data-value", value),
            );
        }
        let zero = self.range.y(0.0, top);
        row.append(
            Line::new()
                .set("data-role", "zero-score-rule")
                .set("x1", PLOT_LEFT)
                .set("x2", PLOT_LEFT + PLOT_WIDTH)
                .set("y1", zero)
                .set("y2", zero)
                .set("stroke", "#81939c")
                .set("stroke-dasharray", "3 3"),
        );
        let mut isolated_points = Group::new();
        for (strand, scores) in [
            (TssStrand::Plus, &self.track.forward_scores),
            (TssStrand::Minus, &self.track.reverse_scores),
        ] {
            isolated_points.append(draw_trace(
                &mut row,
                scores,
                self.range,
                axis,
                top,
                &self.color,
                strand,
                geometry.strand,
                clip,
            ));
        }
        let tss = axis.x(geometry.upstream_bp as f64);
        row.append(
            Line::new()
                .set("data-role", "tss-rule-halo")
                .set("x1", tss)
                .set("x2", tss)
                .set("y1", top)
                .set("y2", top + PLOT_HEIGHT)
                .set("stroke", "#ffffff")
                .set("stroke-width", 5.5),
        );
        row.append(
            Line::new()
                .set("data-role", "tss-rule")
                .set("data-relative-bp", 0)
                .set("x1", tss)
                .set("x2", tss)
                .set("y1", top)
                .set("y2", top + PLOT_HEIGHT)
                .set("stroke", "#172e3c")
                .set("stroke-width", 2.5)
                .add(Title::new("0 bp: annotated transcript-start candidate")),
        );
        row.append(isolated_points);
        self.notes.draw(
            &mut row,
            PLOT_LEFT,
            top + PLOT_HEIGHT + 24.0,
            "matrix-score-notes",
        );
        row
    }
}

fn peak_label(peak: Option<&TssPeak>, geometry: &TssGeometry) -> String {
    peak.map(|peak| {
        format!(
            "{} at {:+} bp",
            number(peak.score),
            geometry
                .relative_at(peak.local_start_0based)
                .expect("validated peak")
        )
    })
    .unwrap_or_else(|| "not reported".into())
}

#[allow(clippy::too_many_arguments)]
fn draw_trace(
    row: &mut Group,
    scores: &[Option<f64>],
    range: ScoreRange,
    axis: LocalAxis,
    top: f64,
    color: &str,
    local_strand: TssStrand,
    genomic_strand: TssStrand,
    clip: bool,
) -> Group {
    let genomic_strand = if local_strand == TssStrand::Plus {
        genomic_strand
    } else {
        genomic_strand.opposite()
    };
    let mut path = String::new();
    let mut connected = false;
    let mut isolated = Vec::new();
    for (index, score) in scores.iter().enumerate() {
        let Some(score) = score else {
            connected = false;
            continue;
        };
        let x = axis.x(index as f64);
        let y = range.y(display_score(*score, clip), top);
        let _ = write!(path, "{} {x:.4} {y:.4} ", if connected { 'L' } else { 'M' });
        if !connected && scores.get(index + 1).is_none_or(Option::is_none) {
            isolated.push((index, x, y));
        }
        connected = true;
    }
    if !path.is_empty() {
        row.append(
            Path::new()
                .set("data-role", "score-trace")
                .set("data-local-strand", local_strand.as_str())
                .set("data-genomic-strand", genomic_strand.as_str())
                .set("data-coordinate", "motif_window_start_0based")
                .set("d", path.trim())
                .set("fill", "none")
                .set("stroke", color)
                .set("stroke-width", 1.8)
                .set("stroke-linejoin", "round")
                .set("stroke-linecap", "round")
                .set(
                    "stroke-dasharray",
                    if local_strand == TssStrand::Minus {
                        "6 4"
                    } else {
                        "none"
                    },
                ),
        );
    }
    let mut markers = Group::new();
    for (index, x, y) in isolated {
        markers.append(
            Circle::new()
                .set("data-role", "isolated-score")
                .set("data-local-start", index)
                .set("data-local-strand", local_strand.as_str())
                .set("data-genomic-strand", genomic_strand.as_str())
                .set("cx", x)
                .set("cy", y)
                .set("r", 2.2)
                .set(
                    "fill",
                    if local_strand == TssStrand::Plus {
                        color
                    } else {
                        "#ffffff"
                    },
                )
                .set("stroke", color)
                .set("stroke-width", 1.5),
        );
    }
    markers
}

fn logo_bits(counts: &[f64; 4]) -> [f64; 4] {
    let max = counts.iter().copied().fold(0.0_f64, f64::max);
    let scaled = counts.map(|count| count / max);
    let total: f64 = scaled.iter().sum();
    let probabilities = scaled.map(|count| count / total);
    let information = (2.0
        + probabilities
            .iter()
            .filter(|p| **p > 0.0)
            .map(|p| p * p.log2())
            .sum::<f64>())
    .clamp(0.0, 2.0);
    probabilities.map(|p| p * information)
}

fn draw_logo(group: &mut Group, matrix: &ResolvedTssMatrix, x: f64, top: f64) {
    // Unit-square outline letters have exact bounds and need no installed font.
    const GLYPHS: [(&str, &str, &str); 4] = [
        (
            "A",
            "#18814b",
            "M0 100 L36 0 H64 L100 100 H76 L68 76 H32 L24 100 Z M40 54 H60 L50 23 Z",
        ),
        (
            "C",
            "#2873bc",
            "M98 16 L82 34 C73 25 66 22 54 22 C34 22 24 32 24 50 C24 68 34 78 54 78 C66 78 74 74 82 65 L98 82 C84 96 72 100 52 100 C18 100 0 81 0 50 C0 19 18 0 52 0 C72 0 86 6 98 16 Z",
        ),
        (
            "G",
            "#a97900",
            "M98 16 L82 34 C73 25 66 22 54 22 C34 22 24 32 24 50 C24 68 34 78 54 78 C63 78 70 76 76 72 V62 H54 V43 H100 V84 C86 96 72 100 52 100 C18 100 0 81 0 50 C0 19 18 0 52 0 C72 0 86 6 98 16 Z",
        ),
        ("T", "#c93f45", "M0 0 H100 V23 H62 V100 H38 V23 H0 Z"),
    ];
    for (segment_index, columns) in matrix.matrix_counts.chunks(LOGO_COLUMNS).enumerate() {
        let y = top + segment_index as f64 * LOGO_SEGMENT_HEIGHT;
        let left = x + 32.0;
        let mut segment = Group::new()
            .set("data-role", "pfm-logo")
            .set("data-accession", matrix.specification.source_id.as_str())
            .set("data-first-position", segment_index * LOGO_COLUMNS + 1)
            .set("data-column-width", LOGO_COLUMN_WIDTH);
        for bits in 0..=2 {
            let yy = y + LOGO_HEIGHT * (1.0 - bits as f64 / 2.0);
            segment.append(
                Line::new()
                    .set("x1", left)
                    .set("x2", left + columns.len() as f64 * LOGO_COLUMN_WIDTH)
                    .set("y1", yy)
                    .set("y2", yy)
                    .set("stroke", "#cad8df")
                    .set("stroke-width", 0.7),
            );
            segment.append(
                text_node(left - 8.0, yy + 4.0, &bits.to_string(), 12.0).set("text-anchor", "end"),
            );
        }
        for (column_index, counts) in columns.iter().enumerate() {
            let position = segment_index * LOGO_COLUMNS + column_index + 1;
            let contributions = logo_bits(counts);
            let mut column = Group::new()
                .set("data-role", "pfm-column")
                .set("data-position-1based", position)
                .set("data-information-bits", contributions.iter().sum::<f64>());
            let mut letters = [0, 1, 2, 3];
            letters.sort_by(|a, b| contributions[*a].total_cmp(&contributions[*b]));
            let mut bottom = y + LOGO_HEIGHT;
            for base in letters {
                let bits = contributions[base];
                if bits <= 0.0 {
                    continue;
                }
                let height = bits * LOGO_HEIGHT / 2.0;
                bottom -= height;
                let (letter, color, path) = GLYPHS[base];
                column.append(
                    Path::new()
                        .set("data-role", "pfm-letter")
                        .set("data-base", letter)
                        .set("data-bits", bits)
                        .set("d", path)
                        .set("fill", color)
                        .set("fill-rule", "evenodd")
                        .set(
                            "transform",
                            format!(
                                "translate({:.4} {:.4}) scale({:.6} {:.9})",
                                left + column_index as f64 * LOGO_COLUMN_WIDTH + 1.0,
                                bottom,
                                (LOGO_COLUMN_WIDTH - 2.0) / 100.0,
                                height / 100.0
                            ),
                        )
                        .add(Title::new(format!(
                            "Position {position}: {letter}, count {}, information {} bits",
                            counts[base],
                            number(bits)
                        ))),
                );
            }
            segment.append(column);
            if column_index == 0 || column_index + 1 == columns.len() || position % 5 == 0 {
                segment.append(
                    text_node(
                        left + (column_index as f64 + 0.5) * LOGO_COLUMN_WIDTH,
                        y + LOGO_HEIGHT + 17.0,
                        &position.to_string(),
                        12.0,
                    )
                    .set("text-anchor", "middle"),
                );
            }
        }
        group.append(segment);
    }
}

fn comparison_block(
    comparison: &TssMatrixComparison,
    geometry: &TssGeometry,
    score_kind: &str,
) -> TextBlock {
    let genomic_strand = if comparison.local_strand == TssStrand::Plus {
        geometry.strand
    } else {
        geometry.strand.opposite()
    };
    let metric = |value: Option<f64>| value.map(number).unwrap_or_else(|| "undefined".into());
    let mut text = format!(
        "{}: {} versus {} | local {} / genomic {}\nPearson: {} | Spearman: {} | paired starts: {} | excluded starts: {}\nInput score: {} | Method: {}",
        comparison.factor_id,
        comparison.left_accession,
        comparison.right_accession,
        comparison.local_strand.as_str(),
        genomic_strand.as_str(),
        metric(comparison.pearson),
        metric(comparison.spearman),
        comparison.paired_window_count,
        comparison.excluded_window_count,
        score_kind,
        if comparison.method.trim().is_empty() {
            "not supplied by report"
        } else {
            &comparison.method
        },
    );
    if comparison.pearson.is_none()
        || comparison.spearman.is_none()
        || comparison.undefined_reason.is_some()
    {
        let _ = write!(
            text,
            "\nUndefined reason: {}",
            comparison
                .undefined_reason
                .as_deref()
                .filter(|s| !s.trim().is_empty())
                .unwrap_or("not supplied by report")
        );
    }
    TextBlock::new(&text, TEXT_WIDTH, 14.0)
}

struct WindowLayout<'a> {
    window: &'a TssProfileWindow,
    intro: TextBlock,
    legend: TextBlock,
    axis: AxisLayout,
    rows: Vec<RowLayout<'a>>,
    comparison_header: TextBlock,
    comparisons: Vec<TextBlock>,
    clip: bool,
}

#[derive(Clone, Debug)]
struct Fragment {
    window: usize,
    rows: Range<usize>,
    comparisons: Range<usize>,
}

impl<'a> WindowLayout<'a> {
    fn new(
        report: &'a TssProfileReport,
        window: &'a TssProfileWindow,
        scale: TssScaleMode,
        position: usize,
        count: usize,
    ) -> Result<Self, String> {
        let record = &window.record;
        let geometry = &record.geometry;
        let evidence = if window.selected {
            window.selection_evidence.as_ref()
        } else {
            None
        };
        let selection_label =
            evidence
                .map(|evidence| evidence.label.as_str())
                .unwrap_or(if window.selected {
                    "Selected"
                } else {
                    "Unselected"
                });
        let mut intro = format!("TSS window {position} of {count} | {selection_label}");
        if let Some(evidence) = evidence {
            let _ = write!(intro, "\n{}", evidence.legend);
        }
        let _ = write!(
            intro,
            "\nPromoter: {}\nAnnotated TSS: {}:{} ({}) | stored interval: {}..{} (ascending)\nTranscripts ({}): {}\nSequence SHA-256: {}",
            record.promoter_id,
            geometry.chromosome,
            geometry.tss_1based,
            geometry.strand.as_str(),
            geometry.start_1based,
            geometry.end_1based,
            record.transcripts.len(),
            if record.transcripts.is_empty() {
                "none supplied".into()
            } else {
                record.transcripts.join(", ")
            },
            record.sequence_sha256,
        );
        let intro = TextBlock::new(&intro, TEXT_WIDTH, 14.0);
        let legend = TextBlock::new(
            &format!(
                "Transcript-oriented sequence: 5'-to-3' left to right; genomic labels {}.\nSolid: motif-local + / genomic {}. Dashed: motif-local - / genomic {}. Black vertical rule: 0 bp, annotated TSS.\nEvery point is a motif-window START on the common sequence axis (not a center or reverse-motif 5' endpoint).\nGray: terminal positions without a full motif window. Amber and path gaps: unavailable windows (including N), not zero scores.",
                if geometry.strand == TssStrand::Minus {
                    "descend"
                } else {
                    "ascend"
                },
                geometry.strand.as_str(),
                geometry.strand.opposite().as_str(),
            ),
            TEXT_WIDTH,
            13.0,
        );
        let panel = &report.panel_resolution.panel;
        let shared = ScoreRange::from_tracks(window.tracks.iter(), panel.clip_negative);
        let by_accession: BTreeMap<_, _> = window
            .tracks
            .iter()
            .map(|track| (track.accession.as_str(), track))
            .collect();
        let rows = report
            .panel_resolution
            .matrices
            .iter()
            .enumerate()
            .map(|(index, matrix)| {
                let track = by_accession[matrix.specification.source_id.as_str()];
                let range = if scale == TssScaleMode::Shared {
                    shared
                } else {
                    ScoreRange::from_tracks(std::iter::once(track), panel.clip_negative)
                };
                RowLayout::new(report, window, matrix, track, range, index)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let comparison_header = TextBlock::new(
            "WITHIN-FACTOR MATRIX COMPARISONS\nStored results on common valid window starts; metric methods retain strand, clipping and smoothing policies. Visual clipping/scaling is not reapplied. Correlation does not identify a biologically correct matrix.",
            TEXT_WIDTH,
            14.0,
        );
        let comparisons = if window.comparisons.is_empty() {
            vec![TextBlock::new(
                "No within-factor comparison pairs reported. A single-matrix factor has zero pairs, not a failed correlation.",
                TEXT_WIDTH,
                14.0,
            )]
        } else {
            window
                .comparisons
                .iter()
                .map(|comparison| comparison_block(comparison, geometry, &panel.score_kind))
                .collect()
        };
        Ok(Self {
            window,
            intro,
            legend,
            axis: AxisLayout::new(geometry),
            rows,
            comparison_header,
            comparisons,
            clip: panel.clip_negative,
        })
    }

    fn header_height(&self) -> f64 {
        self.intro.height() + 12.0 + 24.0 + self.legend.height() + 16.0 + self.axis.height()
    }

    fn fragment_height(&self, fragment: &Fragment) -> f64 {
        self.header_height()
            + self.rows[fragment.rows.clone()]
                .iter()
                .map(|row| row.height + ROW_GAP)
                .sum::<f64>()
            + if fragment.comparisons.is_empty() {
                0.0
            } else {
                self.comparison_header.height()
                    + 12.0
                    + self.comparisons[fragment.comparisons.clone()]
                        .iter()
                        .map(|block| block.height() + 16.0)
                        .sum::<f64>()
            }
    }

    fn fragments(&self, window: usize, capacity: f64) -> Result<Vec<Fragment>, String> {
        let mut output = Vec::new();
        let mut fragment = Fragment {
            window,
            rows: 0..0,
            comparisons: 0..0,
        };
        for row in 0..self.rows.len() {
            let mut candidate = fragment.clone();
            candidate.rows.end = row + 1;
            if self.fragment_height(&candidate) > capacity && !fragment.rows.is_empty() {
                output.push(fragment);
                fragment = Fragment {
                    window,
                    rows: row..row,
                    comparisons: 0..0,
                };
                candidate = fragment.clone();
                candidate.rows.end = row + 1;
            }
            if self.fragment_height(&candidate) > capacity {
                return Err(format!(
                    "TSS {} matrix {} cannot fit a readable row within the 10000-pixel page bound",
                    self.window.record.promoter_id, self.rows[row].matrix.specification.source_id
                ));
            }
            fragment = candidate;
        }
        for comparison in 0..self.comparisons.len() {
            let mut candidate = fragment.clone();
            candidate.comparisons.end = comparison + 1;
            if self.fragment_height(&candidate) > capacity
                && (!fragment.rows.is_empty() || !fragment.comparisons.is_empty())
            {
                output.push(fragment);
                fragment = Fragment {
                    window,
                    rows: self.rows.len()..self.rows.len(),
                    comparisons: comparison..comparison,
                };
                candidate = fragment.clone();
                candidate.comparisons.end = comparison + 1;
            }
            if self.fragment_height(&candidate) > capacity {
                return Err(format!(
                    "TSS {} comparison/header cannot fit the 10000-pixel page bound",
                    self.window.record.promoter_id
                ));
            }
            fragment = candidate;
        }
        output.push(fragment);
        Ok(output)
    }

    fn draw(&self, fragment: &Fragment, top: f64) -> Group {
        let mut group = Group::new()
            .set("data-role", "tss-panel")
            .set("data-promoter-id", self.window.record.promoter_id.as_str())
            .set("data-selected", self.window.selected.to_string())
            .set("data-y", top)
            .set("data-height", self.fragment_height(fragment));
        self.intro.draw(&mut group, MARGIN, top, "tss-identity");
        let mut y = top + self.intro.height() + 12.0;
        let caption = if fragment.rows.is_empty() {
            "Comparison continuation; accession identities repeated below".into()
        } else {
            format!(
                "Matrix rows {}..{} of {}{}",
                fragment.rows.start + 1,
                fragment.rows.end,
                self.rows.len(),
                if fragment.rows.start > 0 || fragment.rows.end < self.rows.len() {
                    " (continued TSS; row sizes unchanged)"
                } else {
                    ""
                }
            )
        };
        group.append(
            text_node(MARGIN, y + 14.0, &caption, 14.0).set("data-role", "track-page-caption"),
        );
        y += 24.0;
        self.legend
            .draw(&mut group, MARGIN, y, "strand-coordinate-legend");
        y += self.legend.height() + 16.0;
        self.axis.draw(&mut group, y);
        y += self.axis.height();
        for row in &self.rows[fragment.rows.clone()] {
            group.append(row.draw(y, self.axis.axis, &self.window.record.geometry, self.clip));
            y += row.height + ROW_GAP;
        }
        if !fragment.comparisons.is_empty() {
            self.comparison_header
                .draw(&mut group, MARGIN, y, "comparison-header");
            y += self.comparison_header.height() + 12.0;
            for comparison in &self.comparisons[fragment.comparisons.clone()] {
                comparison.draw(&mut group, MARGIN, y, "matrix-comparison");
                y += comparison.height() + 16.0;
            }
        }
        group
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gentle_protocol::tss_profiles::{
        JasparPanelTrack, JasparTargetPanel, PANEL_SCHEMA, TssPanelResolution, TssRecord,
        TssReference, TssSelectionEvidence, TssStrandPolicy,
    };
    use svg::node::Attributes;
    use svg::parser::Event;

    // Hand-crafted synthetic reports only, recreated by fixture(). No source
    // sequences, biological observations, real accession PFMs or external files
    // are used. They exercise this renderer, not the engine's scoring formulas.
    fn fixture(upstream: usize, downstream: usize, matrices: usize) -> TssProfileReport {
        let length = upstream + 1 + downstream;
        let matrices = (0..matrices)
            .map(|index| {
                let counts = [[8.0, 0.0, 0.0, 0.0], [6.0, 2.0, 1.0, 1.0], [1.0; 4]];
                ResolvedTssMatrix {
                    specification: JasparPanelTrack {
                        source_id: format!("MA{:04}.1", 9000 + index),
                        factor_id: "SYNTHETIC_FACTOR".into(),
                        label: format!("Synthetic matrix {}", index + 1),
                        display_order: index + 1,
                        color_hint: Some(COLORS[index % COLORS.len()].into()),
                        score_kind: None,
                        track_id: None,
                        provider_kind: None,
                        factor_label: None,
                    },
                    version: "synthetic-v1".into(),
                    consensus: "NEVER_RENDER_THIS_CONSENSUS".into(),
                    matrix_counts: (0..3 + index % 3)
                        .map(|column| counts[column % 3])
                        .collect(),
                    matrix_sha256: "b".repeat(64),
                }
            })
            .collect::<Vec<_>>();
        let tracks = matrices.iter().map(|matrix| {
            let count = length.checked_sub(matrix.matrix_counts.len()).map_or(0, |n| n + 1);
            TssProfileTrack {
                accession: matrix.specification.source_id.clone(), motif_length_bp: matrix.matrix_counts.len(),
                forward_scores: (0..count).map(|i| Some((i % 5) as f64 - 1.0)).collect(),
                reverse_scores: (0..count).map(|i| Some((i % 3) as f64 - 0.5)).collect(),
                forward_maximum: None, reverse_maximum: None, forward_peaks: vec![], reverse_peaks: vec![],
                normalization_reference: serde_json::json!({"background_model":"synthetic IID", "pseudocount":0.1, "quantization_bits":0.1}),
            }
        }).collect();
        TssProfileReport {
            schema: REPORT_SCHEMA.into(),
            reference: TssReference {
                genome_id: "synthetic_genome".into(),
                assembly: "synthetic_assembly".into(),
                annotation_release: None,
            },
            panel_resolution: TssPanelResolution {
                panel: JasparTargetPanel {
                    schema: PANEL_SCHEMA.into(),
                    panel_id: "synthetic_panel".into(),
                    label: "Synthetic renderer tests".into(),
                    score_kind: "llr_bits".into(),
                    clip_negative: false,
                    scale_mode: TssScaleMode::Independent,
                    strand_policy: TssStrandPolicy::Both,
                    calibration_state: TssCalibrationState::MatrixSpecific,
                    calibration_statement:
                        "Synthetic matrix-specific scores; no cross-matrix calibration.".into(),
                    calibration_id: None,
                    calibration_sha256: None,
                    top_hit_count: 10,
                    factors: matrices
                        .iter()
                        .map(|matrix| matrix.specification.clone())
                        .collect(),
                },
                panel_sha256: "a".repeat(64),
                registry_sources: vec![],
                registry_source_url: None,
                matrices,
            },
            inputs: vec![],
            windows: vec![TssProfileWindow {
                record: TssRecord {
                    promoter_id: "synthetic_tss_1".into(),
                    gene_id: "SYNTH_GENE_ID".into(),
                    gene_symbol: "SYNTH_GENE".into(),
                    geometry: TssGeometry {
                        chromosome: "chrSynthetic".into(),
                        strand: TssStrand::Plus,
                        tss_1based: 10001,
                        start_1based: 10001 - upstream as u64,
                        end_1based: 10001 + downstream as u64,
                        upstream_bp: upstream,
                        downstream_bp: downstream,
                    },
                    transcripts: vec![
                        "synthetic_transcript_1".into(),
                        "synthetic_transcript_2".into(),
                    ],
                    sequence_sha256: "c".repeat(64),
                },
                selected: true,
                selection_evidence: None,
                tracks,
                comparisons: vec![],
            }],
            source: None,
            producer_revision: "synthetic".into(),
            producer_executable_sha256: None,
            lockfile_sha256: "d".repeat(64),
            score_policy: BTreeMap::from([(
                "clipping_for_comparison".into(),
                "none (raw scores)".into(),
            )]),
            verification: "Synthetic geometry only; not reference-authenticated".into(),
            warnings: vec![],
            non_claims: NON_CLAIMS.into(),
        }
    }

    fn reverse_geometry(report: &mut TssProfileReport) {
        for window in &mut report.windows {
            let geometry = &mut window.record.geometry;
            geometry.strand = TssStrand::Minus;
            geometry.start_1based = geometry.tss_1based - geometry.downstream_bp as u64;
            geometry.end_1based = geometry.tss_1based + geometry.upstream_bp as u64;
        }
    }

    fn render(report: &TssProfileReport) -> Vec<TssRenderedPage> {
        render_tss_profile_pages(report, &TssProfileRenderOptions::default()).unwrap()
    }

    fn tags(svg: &str, role: &str) -> Vec<Attributes> {
        svg::read(svg)
            .unwrap()
            .filter_map(|event| match event {
                Event::Error(error) => panic!("invalid SVG: {error}"),
                Event::Tag(name, _, attributes)
                    if (role == "svg" && name == "svg")
                        || attributes
                            .get("data-role")
                            .is_some_and(|value| value.to_string() == role) =>
                {
                    Some(attributes)
                }
                _ => None,
            })
            .collect()
    }

    fn numeric(attributes: &Attributes, name: &str) -> f64 {
        attributes[name].to_string().parse().unwrap()
    }

    fn path_points(attributes: &Attributes) -> Vec<(char, f64, f64)> {
        let data = attributes["d"].to_string();
        let tokens: Vec<_> = data.split_whitespace().collect();
        tokens
            .chunks_exact(3)
            .map(|parts| {
                (
                    parts[0].chars().next().unwrap(),
                    parts[1].parse().unwrap(),
                    parts[2].parse().unwrap(),
                )
            })
            .collect()
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!((actual - expected).abs() < 0.0002, "{actual} != {expected}");
    }

    #[test]
    fn plus_minus_geometry_and_nondefault_axes_use_protocol_coordinates() {
        for (upstream, downstream) in [(500, 200), (7, 2), (1, 999), (999, 1), (0, 0)] {
            let mut report = fixture(upstream, downstream, 2);
            for strand in [TssStrand::Plus, TssStrand::Minus] {
                if strand == TssStrand::Minus {
                    reverse_geometry(&mut report);
                }
                let page = &render(&report)[0];
                let root = &tags(&page.svg, "svg")[0];
                assert_close(numeric(root, "width"), 1400.0);
                assert_close(numeric(root, "data-gentle-plot-left"), 255.0);
                assert_close(numeric(root, "data-gentle-plot-right"), 1050.0);
                let geometry = &report.windows[0].record.geometry;
                let axis = LocalAxis {
                    length: upstream + 1 + downstream,
                };
                let ticks = tags(&page.svg, "genomic-tick");
                for index in [0, upstream, upstream + downstream] {
                    let tick = ticks
                        .iter()
                        .find(|tick| numeric(tick, "data-local-start") == index as f64)
                        .unwrap();
                    assert_eq!(
                        tick["data-value"].to_string(),
                        geometry.genomic_at(index).unwrap().to_string()
                    );
                    assert_close(numeric(tick, "data-axis-x"), axis.x(index as f64));
                }
                let mut coordinates: Vec<_> = ticks
                    .iter()
                    .map(|tick| numeric(tick, "data-value"))
                    .collect();
                coordinates.dedup();
                assert!(
                    coordinates
                        .windows(2)
                        .all(|pair| if strand == TssStrand::Plus {
                            pair[0] < pair[1]
                        } else {
                            pair[0] > pair[1]
                        })
                );
                let rules = tags(&page.svg, "tss-rule");
                assert_eq!(rules.len(), 2);
                for rule in rules {
                    assert_close(numeric(&rule, "x1"), axis.x(upstream as f64));
                    assert_close(numeric(&rule, "x2"), axis.x(upstream as f64));
                    assert_close(numeric(&rule, "y2") - numeric(&rule, "y1"), PLOT_HEIGHT);
                    assert!(numeric(&rule, "stroke-width") >= 2.0);
                }
                for trace in tags(&page.svg, "score-trace") {
                    let genomic = if trace["data-local-strand"].to_string() == "+" {
                        strand
                    } else {
                        strand.opposite()
                    };
                    assert_eq!(trace["data-genomic-strand"].to_string(), genomic.as_str());
                }
            }
        }
    }

    #[test]
    fn last_motif_window_does_not_stretch_to_sequence_end() {
        let report = fixture(5, 3, 3);
        let page = &render(&report)[0];
        let traces = tags(&page.svg, "score-trace");
        let terminal = tags(&page.svg, "terminal-unscored");
        assert!(
            page.svg.rfind("data-role=\"tss-rule\"").unwrap()
                > page.svg.rfind("data-role=\"score-trace\"").unwrap()
        );
        for (index, track) in report.windows[0].tracks.iter().enumerate() {
            let last_start = 9 - track.motif_length_bp;
            for trace in &traces[index * 2..index * 2 + 2] {
                let points = path_points(trace);
                assert_close(points[0].1, PLOT_LEFT);
                assert_close(
                    points.last().unwrap().1,
                    PLOT_LEFT + last_start as f64 / 8.0 * PLOT_WIDTH,
                );
                assert!(points.last().unwrap().1 < PLOT_LEFT + PLOT_WIDTH);
                let scored_anchor = 5.min(last_start);
                assert_close(
                    points[scored_anchor].1,
                    PLOT_LEFT + scored_anchor as f64 / 8.0 * PLOT_WIDTH,
                );
            }
            assert_eq!(
                numeric(&terminal[index], "data-local-start"),
                (last_start + 1) as f64
            );
            assert_eq!(numeric(&terminal[index], "data-local-end-exclusive"), 9.0);
            assert!(
                numeric(&terminal[index], "x") > PLOT_LEFT + last_start as f64 / 8.0 * PLOT_WIDTH
            );
            assert_close(
                numeric(&terminal[index], "x") + numeric(&terminal[index], "width"),
                PLOT_LEFT + PLOT_WIDTH,
            );
        }
    }

    #[test]
    fn independent_scales_show_a_hundredfold_difference_without_changing_scores() {
        let mut report = fixture(5, 3, 2);
        for (index, track) in report.windows[0].tracks.iter_mut().enumerate() {
            for score in track
                .forward_scores
                .iter_mut()
                .chain(&mut track.reverse_scores)
            {
                *score = Some(if index == 0 { 2.0 } else { 200.0 });
            }
        }
        let before = serde_json::to_value(&report).unwrap();
        let page = &render(&report)[0];
        let rows = tags(&page.svg, "matrix-row");
        assert_eq!(numeric(&rows[0], "data-y-min"), 0.0);
        assert_eq!(numeric(&rows[0], "data-y-max"), 2.0);
        assert_eq!(numeric(&rows[1], "data-y-max"), 200.0);
        let ticks = tags(&page.svg, "score-tick");
        assert_eq!(
            ticks
                .iter()
                .map(|tick| numeric(tick, "data-value"))
                .collect::<Vec<_>>(),
            vec![0.0, 1.0, 2.0, 0.0, 100.0, 200.0]
        );
        assert_eq!(before, serde_json::to_value(&report).unwrap());
        assert!(page.svg.contains("Equal trace heights"));
    }

    #[test]
    fn shared_ranges_require_typed_exact_calibration_not_units_or_prose() {
        let mut report = fixture(5, 3, 2);
        let shared = TssProfileRenderOptions {
            scale_mode: Some(TssScaleMode::Shared),
            panels_per_page: 1,
        };
        assert!(
            render_tss_profile_pages(&report, &shared)
                .unwrap_err()
                .contains("CrossSourceCalibrated")
        );
        report.panel_resolution.panel.calibration_id = Some("synthetic-calibration-v1".into());
        report.panel_resolution.panel.calibration_sha256 = Some("a".repeat(64));
        assert!(render_tss_profile_pages(&report, &shared).is_err());
        report.panel_resolution.panel.calibration_state =
            TssCalibrationState::CrossSourceCalibrated;
        let page = &render_tss_profile_pages(&report, &shared).unwrap()[0];
        let rows = tags(&page.svg, "matrix-row");
        assert_eq!(rows[0]["data-y-min"], rows[1]["data-y-min"]);
        assert_eq!(rows[0]["data-y-max"], rows[1]["data-y-max"]);
        assert!(page.svg.contains("synthetic-calibration-v1"));
        assert!(page.svg.contains(&"a".repeat(64)));
        for digest in [
            "A".repeat(64),
            format!("sha256:{}", "a".repeat(64)),
            "g".repeat(64),
            "".into(),
        ] {
            let mut invalid = report.clone();
            invalid.panel_resolution.panel.calibration_sha256 = Some(digest);
            assert!(render_tss_profile_pages(&invalid, &shared).is_err());
        }
        for id in [None, Some("".into()), Some(" untrimmed ".into())] {
            let mut invalid = report.clone();
            invalid.panel_resolution.panel.calibration_id = id;
            assert!(render_tss_profile_pages(&invalid, &shared).is_err());
        }
        report.panel_resolution.panel.calibration_statement = "  ".into();
        assert!(render_tss_profile_pages(&report, &shared).is_err());
        let mut panel_default = fixture(5, 3, 2);
        panel_default.panel_resolution.panel.scale_mode = TssScaleMode::Shared;
        assert!(
            render_tss_profile_pages(&panel_default, &TssProfileRenderOptions::default()).is_err()
        );
        assert!(
            render_tss_profile_pages(
                &panel_default,
                &TssProfileRenderOptions {
                    scale_mode: Some(TssScaleMode::Independent),
                    panels_per_page: 1
                }
            )
            .is_ok()
        );
    }

    #[test]
    fn negative_zero_and_null_scores_have_distinct_geometry() {
        let mut report = fixture(5, 3, 1);
        let scores = vec![
            Some(-2.0),
            None,
            Some(0.0),
            Some(4.0),
            None,
            Some(-1.0),
            Some(0.0),
        ];
        report.windows[0].tracks[0].forward_scores = scores.clone();
        report.windows[0].tracks[0].reverse_scores = vec![None; 7];
        let raw = render(&report).remove(0);
        assert_eq!(
            numeric(&tags(&raw.svg, "matrix-row")[0], "data-y-min"),
            -2.0
        );
        let trace = tags(&raw.svg, "score-trace").remove(0);
        assert_eq!(
            path_points(&trace)
                .iter()
                .filter(|point| point.0 == 'M')
                .count(),
            3
        );
        assert_eq!(tags(&raw.svg, "isolated-score").len(), 1);
        assert_eq!(tags(&raw.svg, "unavailable-windows").len(), 3);
        report.panel_resolution.panel.clip_negative = true;
        let clipped = render(&report).remove(0);
        let row = tags(&clipped.svg, "matrix-row").remove(0);
        assert_eq!(numeric(&row, "data-y-min"), 0.0);
        let trace = tags(&clipped.svg, "score-trace").remove(0);
        assert_close(
            path_points(&trace)[0].2,
            numeric(&row, "data-y") + PLOT_TOP + PLOT_HEIGHT,
        );
        assert_eq!(report.windows[0].tracks[0].forward_scores, scores);
        assert!(!clipped.svg.contains("NaN"));
    }

    #[test]
    fn zero_negative_missing_extreme_and_too_long_motifs_remain_finite() {
        for value in [
            Some(0.0),
            Some(-3.0),
            None,
            Some(f64::MAX),
            Some(-f64::MAX),
            Some(1e-12),
        ] {
            let mut report = fixture(5, 3, 1);
            report.windows[0].tracks[0].forward_scores.fill(value);
            report.windows[0].tracks[0].reverse_scores.fill(value);
            let page = &render(&report)[0];
            let row = tags(&page.svg, "matrix-row").remove(0);
            assert!(numeric(&row, "data-y-min").is_finite());
            assert!(numeric(&row, "data-y-max").is_finite());
            for trace in tags(&page.svg, "score-trace") {
                assert!(
                    path_points(&trace)
                        .iter()
                        .all(|point| point.1.is_finite() && point.2.is_finite())
                );
            }
            if value.is_none() {
                assert!(tags(&page.svg, "score-trace").is_empty());
                assert!(page.svg.contains("No valid scored windows"));
            }
        }
        let report = fixture(0, 0, 1);
        let page = &render(&report)[0];
        let terminal = tags(&page.svg, "terminal-unscored").remove(0);
        assert_eq!(numeric(&terminal, "width"), PLOT_WIDTH);
        assert!(tags(&page.svg, "score-trace").is_empty());
        let mut single = fixture(0, 0, 1);
        single.panel_resolution.matrices[0]
            .matrix_counts
            .truncate(1);
        let track = &mut single.windows[0].tracks[0];
        track.motif_length_bp = 1;
        track.forward_scores = vec![Some(0.0)];
        track.reverse_scores = vec![None];
        let svg = &render(&single)[0].svg;
        assert_eq!(tags(svg, "isolated-score").len(), 1);
        assert!(tags(svg, "terminal-unscored").is_empty());
        assert!(
            svg.rfind("data-role=\"isolated-score\"").unwrap()
                > svg.rfind("data-role=\"tss-rule\"").unwrap()
        );
    }

    #[test]
    fn pfm_logos_are_information_content_stacks_never_consensus() {
        assert_eq!(logo_bits(&[10.0, 0.0, 0.0, 0.0]), [2.0, 0.0, 0.0, 0.0]);
        assert_eq!(logo_bits(&[1.0; 4]), [0.0; 4]);
        let mixed = logo_bits(&[6.0, 2.0, 1.0, 1.0]);
        assert_close(mixed.iter().sum(), 0.4290494055);
        assert_close(mixed[0], 0.2574296433);
        assert_eq!(logo_bits(&[f64::MAX; 4]), [0.0; 4]);
        let mut report = fixture(80, 20, 1);
        report.panel_resolution.matrices[0].matrix_counts = vec![[6.0, 2.0, 1.0, 1.0]; 45];
        report.windows[0].tracks[0].motif_length_bp = 45;
        report.windows[0].tracks[0].forward_scores = vec![Some(0.0); 57];
        report.windows[0].tracks[0].reverse_scores = vec![Some(0.0); 57];
        let page = &render(&report)[0];
        let logos = tags(&page.svg, "pfm-logo");
        assert_eq!(logos.len(), 3);
        assert_eq!(
            logos
                .iter()
                .map(|logo| numeric(logo, "data-first-position"))
                .collect::<Vec<_>>(),
            vec![1.0, 21.0, 41.0]
        );
        assert!(
            logos
                .iter()
                .all(|logo| numeric(logo, "data-column-width") >= 8.0)
        );
        assert_eq!(tags(&page.svg, "pfm-letter").len(), 45 * 4);
        assert!(!page.svg.contains("NEVER_RENDER_THIS_CONSENSUS"));
    }

    #[test]
    fn resolution_supplies_row_order_labels_colors_and_exact_score_units() {
        let mut report = fixture(5, 3, 3);
        report.windows[0].tracks.reverse();
        report.panel_resolution.matrices.swap(0, 2);
        report.panel_resolution.matrices[0].specification.label =
            "Resolution label, not a track name".into();
        report.panel_resolution.matrices[0].specification.color_hint = Some("#123456".into());
        for kind in [
            "llr_bits",
            "true_log_odds_bits",
            "llr_quantile",
            "true_log_odds_quantile",
            "llr_background_quantile",
            "true_log_odds_background_quantile",
            "llr_background_tail_log10",
            "true_log_odds_background_tail_log10",
        ] {
            report.panel_resolution.panel.score_kind = kind.into();
            let page = &render(&report)[0];
            let rows = tags(&page.svg, "matrix-row");
            assert_eq!(rows[0]["data-accession"].to_string(), "MA9002.1");
            assert_eq!(rows[1]["data-accession"].to_string(), "MA9001.1");
            assert_eq!(
                tags(&page.svg, "score-trace")[0]["stroke"].to_string(),
                "#123456"
            );
            if kind.ends_with("tail_log10") {
                assert!(page.svg.contains("-log10(modeled"));
                assert!(!page.svg.contains("bits (true log odds)"));
            }
        }
    }

    #[test]
    fn long_labels_and_near_endpoint_tss_have_no_gutter_or_tick_collisions() {
        let mut report = fixture(1, 999, 2);
        report.panel_resolution.matrices[0].specification.label = format!(
            "Long <escaped> & label {} {}",
            "W".repeat(280),
            "\u{00c4}".repeat(60)
        );
        report.panel_resolution.matrices[0].specification.factor_id = "FACTOR".repeat(30);
        report.windows[0]
            .record
            .transcripts
            .push("synthetic_transcript_".repeat(50));
        let layout =
            WindowLayout::new(&report, &report.windows[0], TssScaleMode::Independent, 1, 1)
                .unwrap();
        assert!(layout.rows[0].label.lines.len() > 5);
        assert!(layout.axis.tiers >= 2);
        for row in &layout.rows {
            for block in [&row.label, &row.identity, &row.notes] {
                assert!(
                    block
                        .lines
                        .iter()
                        .all(|line| text_width(line, block.size) <= block.width + 0.0001)
                );
            }
            assert!(16.0 + row.label.height() < row.logo_top);
            assert!(row.logo_top + LOGO_SEGMENT_HEIGHT <= row.height);
            assert!(PLOT_TOP + PLOT_HEIGHT + 24.0 + row.notes.height() <= row.height);
        }
        for tier in 0..layout.axis.tiers {
            let ticks: Vec<_> = layout
                .axis
                .ticks
                .iter()
                .filter(|tick| tick.tier == tier)
                .collect();
            for pair in ticks.windows(2) {
                assert!(
                    pair[0].center + pair[0].width * 0.5 + 15.0
                        < pair[1].center - pair[1].width * 0.5
                );
            }
        }
        let page = &render(&report)[0];
        assert!(page.svg.contains("&lt;escaped&gt;"));
        assert!(page.svg.contains("&amp;"));
        let height = numeric(&tags(&page.svg, "svg")[0], "height");
        for event in svg::read(&page.svg).unwrap() {
            if let Event::Tag("text", _, attributes) = event
                && attributes.contains_key("x")
            {
                let width = numeric(&attributes, "textLength");
                let x = numeric(&attributes, "x");
                let left = match attributes
                    .get("text-anchor")
                    .map(|value| value.to_string())
                    .as_deref()
                {
                    Some("middle") => x - width * 0.5,
                    Some("end") => x - width,
                    _ => x,
                };
                assert!(
                    left >= 0.0 && left + width <= PAGE_WIDTH,
                    "text outside page: {attributes:?}"
                );
                assert!(
                    numeric(&attributes, "y") + numeric(&attributes, "font-size") * 0.3 < height
                );
                assert!(numeric(&attributes, "font-size") >= 12.0);
            }
        }
        for tick in tags(&page.svg, "score-tick") {
            assert!(numeric(&tick, "x") >= PLOT_RIGHT);
            assert!(numeric(&tick, "x") + numeric(&tick, "textLength") <= PAGE_WIDTH - MARGIN);
        }
    }

    #[test]
    fn thirty_rows_fit_a_readable_default_page() {
        let report = fixture(500, 200, 30);
        let pages = render(&report);
        assert_eq!(pages.len(), 1);
        let rows = tags(&pages[0].svg, "matrix-row");
        assert_eq!(rows.len(), 30);
        assert!(
            rows.iter()
                .all(|row| numeric(row, "data-height") >= MIN_ROW_HEIGHT)
        );
        assert!(numeric(&tags(&pages[0].svg, "svg")[0], "height") <= MAX_PAGE_HEIGHT);
    }

    #[test]
    fn thirty_one_tss_pages_group_without_sorting_or_dropping_selected_records() {
        let mut report = fixture(5, 3, 1);
        let original = report.windows[0].clone();
        report.windows = (0..31)
            .map(|i| {
                let mut window = original.clone();
                window.record.promoter_id = format!("synthetic_tss_{i:02}");
                window.selected = i % 2 == 1;
                window
            })
            .collect();
        let pages = render(&report);
        assert_eq!(pages.len(), 31);
        for (index, page) in pages.iter().enumerate() {
            assert_eq!(page.page_number, index + 1);
            assert_eq!(page.page_count, 31);
            assert_eq!(page.promoter_ids, vec![format!("synthetic_tss_{index:02}")]);
        }
        let grouped = render_tss_profile_pages(
            &report,
            &TssProfileRenderOptions {
                scale_mode: None,
                panels_per_page: 3,
            },
        )
        .unwrap();
        assert_eq!(grouped.len(), 11);
        assert_eq!(grouped[10].promoter_ids, vec!["synthetic_tss_30"]);
        let bounded = render_tss_profile_pages(
            &report,
            &TssProfileRenderOptions {
                scale_mode: None,
                panels_per_page: 32,
            },
        )
        .unwrap();
        assert!(bounded.len() > 1 && bounded.len() < 31);
        assert_eq!(
            bounded
                .iter()
                .flat_map(|page| &page.promoter_ids)
                .collect::<Vec<_>>(),
            pages
                .iter()
                .flat_map(|page| &page.promoter_ids)
                .collect::<Vec<_>>()
        );
        assert!(
            bounded
                .iter()
                .all(|page| numeric(&tags(&page.svg, "svg")[0], "height") <= MAX_PAGE_HEIGHT)
        );
    }

    fn synthetic_selection_evidence(factor: &str) -> TssSelectionEvidence {
        TssSelectionEvidence {
            label: format!(
                "Selected in integrated report \u{2014} {factor} CUT&RUN-supported TSS window"
            ),
            legend: format!(
                "The TSS belonged to an integrated report panel under a descriptive {factor} CUT&RUN window criterion. Selection does not establish TSS usage, direct binding, or activity."
            ),
            criterion: "synthetic_descriptive_window_membership".into(),
            factor: Some(factor.into()),
        }
    }

    fn rendered_text(svg: &str) -> String {
        let text = svg::read(svg)
            .unwrap()
            .filter_map(|event| match event {
                Event::Text(text) => Some(text),
                Event::Error(error) => panic!("invalid SVG: {error}"),
                _ => None,
            })
            .collect::<Vec<_>>()
            .join(" ");
        text.split_whitespace().collect::<Vec<_>>().join(" ")
    }

    fn assert_selection_evidence_count(svg: &str, evidence: &TssSelectionEvidence, count: usize) {
        let text = rendered_text(svg);
        for value in [&evidence.label, &evidence.legend] {
            let escaped = svg::node::Text::new(value.as_str()).to_string();
            assert_eq!(
                text.matches(escaped.as_str()).count(),
                count,
                "selection text: {value}"
            );
        }
    }

    #[test]
    fn missing_annotation_release_is_not_inferred_from_the_genome_label() {
        let mut report = fixture(5, 3, 1);
        report.reference.genome_id = "synthetic_genome_release_999".into();
        assert!(report.reference.annotation_release.is_none());
        assert!(report.source.is_none());
        assert!(report.producer_executable_sha256.is_none());
        let page = &render(&report)[0];
        let text = rendered_text(&page.svg);
        assert!(text.contains("annotation: not separately declared"));
        assert!(!text.contains("annotation: 999"));
        report.reference.annotation_release = Some("explicit_synthetic_annotation_123".into());
        let text = rendered_text(&render(&report)[0].svg);
        assert!(text.contains("annotation: explicit_synthetic_annotation_123"));
        assert!(!text.contains("annotation: not separately declared"));
    }

    #[test]
    fn selection_evidence_repeats_on_every_track_continuation_page() {
        let mut report = fixture(5, 3, 90);
        reverse_geometry(&mut report);
        let evidence = synthetic_selection_evidence("SYNTH_SELECTION_A");
        report.windows[0].selection_evidence = Some(evidence.clone());
        let before = serde_json::to_value(&report).unwrap();
        let pages = render(&report);
        assert!(pages.len() > 1);
        for page in &pages {
            assert_selection_evidence_count(&page.svg, &evidence, 1);
            assert!(rendered_text(&page.svg).contains("annotation: not separately declared"));
            assert!(numeric(&tags(&page.svg, "svg")[0], "height") <= MAX_PAGE_HEIGHT);
            assert_eq!(page.promoter_ids, vec!["synthetic_tss_1"]);
        }
        assert_eq!(before, serde_json::to_value(&report).unwrap());
    }

    #[test]
    fn selection_evidence_repeats_on_comparison_only_pages() {
        let mut report = fixture(5, 3, 2);
        let evidence = synthetic_selection_evidence("SYNTH_SELECTION_B");
        report.windows[0].selection_evidence = Some(evidence.clone());
        // Repeated synthetic rows stress pagination, not correlation computation.
        report.windows[0].comparisons = vec![
            TssMatrixComparison {
                factor_id: "SYNTHETIC_FACTOR".into(),
                left_accession: "MA9000.1".into(),
                right_accession: "MA9001.1".into(),
                local_strand: TssStrand::Plus,
                paired_window_count: 5,
                excluded_window_count: 2,
                pearson: None,
                spearman: None,
                undefined_reason: Some("constant_signal".into()),
                method: "unclipped raw scores; no smoothing; forward-forward".into(),
            };
            140
        ];
        let pages = render(&report);
        assert!(
            pages
                .iter()
                .any(|page| tags(&page.svg, "matrix-row").is_empty())
        );
        for page in &pages {
            assert_selection_evidence_count(&page.svg, &evidence, 1);
            assert!(numeric(&tags(&page.svg, "svg")[0], "height") <= MAX_PAGE_HEIGHT);
        }
    }

    #[test]
    fn grouped_selection_labels_are_data_driven_and_do_not_leak_to_unselected_panels() {
        let mut report = fixture(5, 3, 1);
        let first_evidence = synthetic_selection_evidence("SYNTH_SELECTION_A");
        let second_evidence = synthetic_selection_evidence("SYNTH_SELECTION_B");
        report.windows[0].selection_evidence = Some(first_evidence.clone());
        let mut second = report.windows[0].clone();
        second.record.promoter_id = "synthetic_tss_2".into();
        second.selection_evidence = Some(second_evidence.clone());
        let mut unselected = report.windows[0].clone();
        unselected.record.promoter_id = "synthetic_tss_3".into();
        unselected.selected = false;
        // Even a stale evidence object cannot label an unselected panel.
        report.windows.extend([second, unselected]);
        let pages = render_tss_profile_pages(
            &report,
            &TssProfileRenderOptions {
                scale_mode: None,
                panels_per_page: 2,
            },
        )
        .unwrap();
        assert_eq!(pages.len(), 2);
        assert_selection_evidence_count(&pages[0].svg, &first_evidence, 1);
        assert_selection_evidence_count(&pages[0].svg, &second_evidence, 1);
        assert_selection_evidence_count(&pages[1].svg, &first_evidence, 0);
        assert_selection_evidence_count(&pages[1].svg, &second_evidence, 0);
        let text = rendered_text(&pages[1].svg);
        assert!(text.contains("| Unselected"));
        assert!(!text.contains("| Selected"));

        report.windows.truncate(1);
        report.windows[0].selection_evidence = None;
        let text = rendered_text(&render(&report)[0].svg);
        assert!(text.contains("| Selected"));
        assert!(!text.contains("Selected in integrated report"));
        assert!(!text.contains("CUT&amp;RUN-supported"));
    }

    #[test]
    fn gene_group_order_and_per_gene_page_counts_follow_first_appearance() {
        let mut report = fixture(5, 3, 1);
        let original = report.windows[0].clone();
        report.windows = ["Z", "A", "Z", "A", "Z"]
            .iter()
            .enumerate()
            .map(|(index, gene)| {
                let mut window = original.clone();
                window.record.gene_id = gene.to_string();
                window.record.gene_symbol = format!("SYMBOL_{gene}");
                window.record.promoter_id = format!("tss_{index}");
                window
            })
            .collect();
        let pages = render(&report);
        assert_eq!(
            pages
                .iter()
                .map(|page| (page.gene_id.as_str(), page.page_number, page.page_count))
                .collect::<Vec<_>>(),
            vec![
                ("Z", 1, 3),
                ("Z", 2, 3),
                ("Z", 3, 3),
                ("A", 1, 2),
                ("A", 2, 2)
            ]
        );
        assert_eq!(
            pages
                .iter()
                .flat_map(|page| page.promoter_ids.iter().map(String::as_str))
                .collect::<Vec<_>>(),
            vec!["tss_0", "tss_2", "tss_4", "tss_1", "tss_3"]
        );
    }

    #[test]
    fn excessive_height_continues_whole_rows_and_repeats_identity_and_footer() {
        let report = fixture(5, 3, 90);
        let pages = render(&report);
        assert!(pages.len() > 1);
        let mut accessions = Vec::new();
        for page in &pages {
            assert!(numeric(&tags(&page.svg, "svg")[0], "height") <= MAX_PAGE_HEIGHT);
            assert_eq!(page.promoter_ids, vec!["synthetic_tss_1"]);
            assert_eq!(tags(&page.svg, "tss-identity").len(), 1);
            assert_eq!(tags(&page.svg, "strand-coordinate-legend").len(), 1);
            assert_eq!(tags(&page.svg, "calibration-nonclaims-footer").len(), 1);
            for row in tags(&page.svg, "matrix-row") {
                assert!(numeric(&row, "data-height") >= MIN_ROW_HEIGHT);
                accessions.push(row["data-accession"].to_string());
            }
        }
        assert_eq!(
            accessions,
            report
                .panel_resolution
                .matrices
                .iter()
                .map(|matrix| matrix.specification.source_id.clone())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn samefactor_comparisons_preserve_methods_and_explicit_undefined_metrics() {
        let mut report = fixture(5, 3, 2);
        reverse_geometry(&mut report);
        report.windows[0].comparisons = vec![TssMatrixComparison {
            factor_id: "SYNTHETIC_FACTOR".into(), left_accession: "MA9000.1".into(), right_accession: "MA9001.1".into(),
            local_strand: TssStrand::Minus, paired_window_count: 5, excluded_window_count: 2,
            pearson: None, spearman: Some(0.0), undefined_reason: Some("constant_signal".into()),
            method: "raw computational scores; no clipping; unsmoothed; paired local starts; reverse-reverse".into(),
        }];
        let before = serde_json::to_value(&report).unwrap();
        let page = &render(&report)[0];
        assert!(page.svg.contains("Pearson: undefined | Spearman: 0"));
        assert!(page.svg.contains("constant_signal"));
        assert!(page.svg.contains("paired starts: 5 | excluded starts: 2"));
        assert!(page.svg.contains("local - / genomic +"));
        assert!(
            page.svg
                .contains("unsmoothed; paired local starts; reverse-reverse")
        );
        assert_eq!(before, serde_json::to_value(&report).unwrap());
        let block = comparison_block(
            &report.windows[0].comparisons[0],
            &report.windows[0].record.geometry,
            "llr_bits",
        )
        .lines;
        report.panel_resolution.panel.clip_negative = true;
        report.panel_resolution.matrices[0].specification.color_hint = Some("#ffffff".into());
        let _ = render(&report);
        assert_eq!(
            block,
            comparison_block(
                &report.windows[0].comparisons[0],
                &report.windows[0].record.geometry,
                "llr_bits"
            )
            .lines
        );
    }

    #[test]
    fn rejects_unsafe_or_inconsistent_reports_without_partial_pages() {
        let report = fixture(5, 3, 2);
        for count in [0, 33, usize::MAX] {
            assert!(
                render_tss_profile_pages(
                    &report,
                    &TssProfileRenderOptions {
                        scale_mode: None,
                        panels_per_page: count
                    }
                )
                .is_err()
            );
        }
        let invalid: Vec<Box<dyn Fn(&mut TssProfileReport)>> = vec![
            Box::new(|r| r.windows[0].tracks[0].forward_scores[0] = Some(f64::NAN)),
            Box::new(|r| {
                r.windows[0].tracks[0].forward_scores.pop();
            }),
            Box::new(|r| {
                r.windows[0].tracks.pop();
            }),
            Box::new(|r| {
                r.windows[0].record.geometry.tss_1based += 1;
            }),
            Box::new(|r| {
                r.panel_resolution.matrices[0].matrix_counts.clear();
            }),
            Box::new(|r| {
                r.panel_resolution.matrices[0].matrix_counts[0] = [0.0; 4];
            }),
            Box::new(|r| {
                r.panel_resolution.matrices[0].matrix_counts[0][0] = -1.0;
            }),
            Box::new(|r| {
                r.panel_resolution.matrices[0].specification.score_kind =
                    Some("true_log_odds_bits".into());
            }),
            Box::new(|r| {
                r.panel_resolution.panel.score_kind = "not_a_score_kind".into();
            }),
            Box::new(|r| {
                r.panel_resolution.matrices[0].specification.color_hint =
                    Some("url(https://example.invalid/color)".into());
            }),
            Box::new(|r| {
                r.panel_resolution.matrices[1].specification.display_order = 1;
            }),
            Box::new(|r| {
                r.windows.push(r.windows[0].clone());
            }),
        ];
        for mutate in invalid {
            let mut bad = report.clone();
            mutate(&mut bad);
            assert!(render_tss_profile_pages(&bad, &TssProfileRenderOptions::default()).is_err());
        }
        let mut empty = report.clone();
        empty.windows.clear();
        assert!(render(&empty).is_empty());
    }

    #[test]
    fn svg_is_byte_deterministic_and_never_omits_mandatory_nonclaims() {
        let mut report = fixture(5, 3, 2);
        report.non_claims.clear();
        assert_eq!(render(&report), render(&report));
        let footer = footer_block(&report, TssScaleMode::Independent)
            .lines
            .join(" ");
        assert!(footer.contains("not measured binding"));
        assert!(footer.contains("Cross-matrix magnitudes are not comparable"));
        assert!(footer.contains("not experimentally established initiation sites"));
        assert!(footer.contains("does not establish autoregulation"));
    }

    #[test]
    #[ignore = "writes synthetic SVGs to the system temp directory for explicit visual QA"]
    fn write_synthetic_visual_examples() {
        let directory = std::env::temp_dir().join("gentle-report-tss-svg-qa");
        std::fs::create_dir_all(&directory).unwrap();
        for (name, mut report) in [
            ("plus", fixture(500, 200, 3)),
            ("minus", fixture(7, 2, 3)),
            ("thirty-rows", fixture(500, 200, 30)),
            ("long-logo", fixture(80, 20, 1)),
        ] {
            if name == "minus" {
                reverse_geometry(&mut report);
            }
            if name == "long-logo" {
                report.panel_resolution.matrices[0].specification.label = "Long label that must wrap instead of colliding with numeric y ticks and the score plot".into();
                report.panel_resolution.matrices[0].matrix_counts = vec![[6.0, 2.0, 1.0, 1.0]; 45];
                report.windows[0].tracks[0].motif_length_bp = 45;
                report.windows[0].tracks[0].forward_scores = vec![Some(0.0); 57];
                report.windows[0].tracks[0].reverse_scores = vec![None; 57];
            }
            for page in render(&report) {
                let path = directory.join(format!("{name}-{}.svg", page.page_number));
                std::fs::write(&path, page.svg).unwrap();
                println!("{}", path.display());
            }
        }
    }
}
