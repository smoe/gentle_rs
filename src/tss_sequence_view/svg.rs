//! Deterministic export of the validated native TSS presentation, without rescoring.

use super::*;
use std::fmt::Write as _;

/// Presentation-only selection. Lane indices refer to the bound view, in its
/// original order; an empty list exports the axis and provenance, not all lanes.
#[derive(Clone, Debug, Serialize)]
pub struct TssViewSvgOptions {
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub lane_indices: Vec<usize>,
    pub width_px: u32,
    pub print_size_mm: Option<(f32, f32)>,
}

fn xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

fn text(svg: &mut String, x: f64, y: f64, s: &str) {
    let _ = write!(
        svg,
        "<text x=\"{x:.2}\" y=\"{y:.2}\" font-size=\"11\">{}</text>",
        xml(s)
    );
}

fn lines(s: &str, width: usize) -> Vec<String> {
    s.chars()
        .collect::<Vec<_>>()
        .chunks(width)
        .map(|c| c.iter().collect())
        .collect()
}

/// Render the selected horizontal span and all enabled lanes, including lanes
/// scrolled vertically out of sight. Source hashes are provenance, not a receipt
/// audit. Refuse excessive work rather than exporting a silently partial plot.
pub fn render_tss_view_svg(
    view: &TssSequenceView,
    options: &TssViewSvgOptions,
) -> Result<String, String> {
    let start = options.start_0based;
    let end = options.end_0based_exclusive;
    let length = view.geometry.length().ok_or("Invalid TSS geometry")?;
    if start >= end
        || end > length
        || !(1000..=5000).contains(&options.width_px)
        || options.lane_indices.len() > 128
        || options.lane_indices.windows(2).any(|w| w[0] >= w[1])
        || options.lane_indices.iter().any(|&i| i >= view.lanes.len())
    {
        return Err("Invalid TSS SVG span, width or ordered lane selection".into());
    }
    let mut work = 0usize;
    for &i in &options.lane_indices {
        let lane = &view.lanes[i];
        if !lane.scale_min.is_finite()
            || !lane.scale_max.is_finite()
            || lane.scale_min >= lane.scale_max
        {
            return Err(format!("Invalid scale for TSS lane {}", lane.id));
        }
        let cells = if let Some(trace) = &lane.trace {
            if trace.forward.len() != trace.reverse.len()
                || trace.forward.len() > length
                || trace
                    .forward
                    .iter()
                    .chain(&trace.reverse)
                    .flatten()
                    .any(|s| !s.is_finite())
            {
                return Err(format!("Invalid score-array geometry for {}", lane.id));
            }
            end.min(trace.forward.len()).saturating_sub(start)
        } else {
            lane.features
                .iter()
                .filter(|f| f.start < end && f.end > start)
                .count()
        };
        work = work.saturating_add(cells);
        if cells > 10_000 || work > 100_000 {
            return Err("TSS SVG paint budget exceeded; zoom in or select fewer lanes".into());
        }
    }
    let width = options.width_px as f64;
    let left = 215.0;
    let right = width - 250.0;
    let x = |p: usize| left + (p as f64 - start as f64) / (end - start) as f64 * (right - left);
    let mut notes = vec![
        "Native TSS presentation; export performs no rescoring or retrieval. All selected lanes, not a screenshot. Unavailable is not zero; stored peaks are not full curves.".into(),
        "Report content/sequence binding is not a full receipt audit. Local model scores and imported raw scores have separate scales.".into(),
        view.provenance.clone(),
    ];
    notes.extend(view.warnings.clone());
    if let Some(local) = &view.local_scoring {
        notes.push(format!("Locally computed {} | report SHA-256 {} | cache {} | producer {}. Model scores, not experimental binding or luciferase activity.", local.request.score_kind.as_str(), local.report_sha256, local.cache_key_sha256, local.producer_revision));
    }
    if let Some(profile) = &view.profile {
        notes.push(format!(
            "Report file SHA-256: {}; producer: {}; panel: {}",
            profile.file_sha256, profile.producer_revision, profile.panel_id
        ));
        notes.extend(profile.warnings.clone());
    }
    let note_lines: Vec<_> = notes
        .iter()
        .flat_map(|s| lines(s, ((width - 40.0) / 7.0) as usize))
        .collect();
    let bottom = 130.0 + options.lane_indices.len() as f64 * 160.0;
    let height = bottom + 30.0 + note_lines.len() as f64 * 15.0;
    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width:.0} {height:.0}\" "
    );
    if let Some((w, h)) = options.print_size_mm {
        if !w.is_finite() || !h.is_finite() || w <= 0.0 || h <= 0.0 {
            return Err("Invalid print dimensions".into());
        }
        let _ = write!(svg, "width=\"{w}mm\" height=\"{h}mm\">");
    } else {
        let _ = write!(svg, "width=\"{width:.0}\" height=\"{height:.0}\">");
    }
    let view_bytes = serde_json::to_vec(view).map_err(|e| e.to_string())?;
    let metadata = serde_json::json!({
        "schema": "gentle.tss_view_svg_projection.v1", "view_sha256": sha256_hex_bytes(&view_bytes),
        "options": options, "sequence_sha256": view.sequence_sha256,
        "geometry": view.geometry, "assembly": view.assembly,
        "genome_id": view.genome_id, "annotation_release": view.annotation_release,
        "profile": view.profile,
        "local_scoring": view.local_scoring,
    });
    let _ = write!(
        svg,
        "<metadata>{}</metadata><rect width=\"100%\" height=\"100%\" fill=\"white\"/><g font-family=\"monospace\" fill=\"#222222\">",
        xml(&metadata.to_string())
    );
    text(
        &mut svg,
        20.0,
        22.0,
        &view
            .title
            .chars()
            .take(((width - 40.0) / 7.0) as usize)
            .collect::<String>(),
    );
    text(
        &mut svg,
        20.0,
        42.0,
        &format!(
            "{} | {} | local {}..{} | {} lanes",
            view.assembly,
            view.geometry.strand.as_str(),
            start + 1,
            end,
            options.lane_indices.len()
        ),
    );
    text(&mut svg, 20.0, 68.0, "TSS-relative bp");
    text(&mut svg, 20.0, 85.0, "Genomic base");
    text(&mut svg, 20.0, 102.0, "Local base");
    let mut ticks = vec![start, start + (end - start - 1) / 2, end - 1];
    ticks.sort_unstable();
    ticks.dedup();
    for p in ticks {
        text(
            &mut svg,
            x(p),
            68.0,
            &format!("{:+}", p as i64 - view.geometry.upstream_bp as i64),
        );
        text(
            &mut svg,
            x(p),
            85.0,
            &view
                .geometry
                .genomic_at(p)
                .ok_or("Invalid genomic coordinate")?
                .to_string(),
        );
        text(&mut svg, x(p), 102.0, &(p + 1).to_string());
    }
    for (row, &index) in options.lane_indices.iter().enumerate() {
        let lane = &view.lanes[index];
        let top = 130.0 + row as f64 * 160.0;
        let low = top + 120.0;
        let y =
            |score: f64| low - (score - lane.scale_min) / (lane.scale_max - lane.scale_min) * 120.0;
        let _ = write!(
            svg,
            "<g data-lane-id=\"{}\"><title>{}</title><defs><clipPath id=\"lane-{row}\"><rect x=\"{left}\" y=\"{top}\" width=\"{}\" height=\"120\"/></clipPath></defs><rect x=\"{left}\" y=\"{top}\" width=\"{}\" height=\"120\" fill=\"#f5f6f7\"/>",
            xml(&lane.id),
            xml(&format!(
                "{}; {}; {}; {}",
                lane.label, lane.units, lane.state, lane.details
            )),
            right - left,
            right - left
        );
        for (n, line) in lines(&lane.label, 25).iter().take(7).enumerate() {
            text(&mut svg, 15.0, top + 14.0 + n as f64 * 14.0, line);
        }
        for (n, line) in lines(&format!("{}; {}", lane.units, lane.state), 32)
            .iter()
            .take(7)
            .enumerate()
        {
            text(&mut svg, right + 15.0, top + 14.0 + n as f64 * 14.0, line);
        }
        text(
            &mut svg,
            left,
            low + 16.0,
            &format!(
                "scale {:.3} .. {:.3}{}",
                lane.scale_min,
                lane.scale_max,
                if lane.trace.as_ref().is_some_and(|t| t.range_is_fallback) {
                    " (display fallback)"
                } else {
                    ""
                }
            ),
        );
        let _ = write!(svg, "<g clip-path=\"url(#lane-{row})\">");
        let tss = view.geometry.upstream_bp;
        if (start..end).contains(&tss) {
            let _ = write!(
                svg,
                "<path data-role=\"tss\" d=\"M {:.2},{top} V {low}\" stroke=\"#555555\"/>",
                x(tss)
            );
        }
        if let Some(trace) = &lane.trace {
            let visible_end = end.min(trace.forward.len());
            if visible_end < end {
                let _ = write!(
                    svg,
                    "<rect data-role=\"terminal-unavailable\" x=\"{:.2}\" y=\"{top}\" width=\"{:.2}\" height=\"120\" fill=\"#dddddd\"><title>No complete motif window starts here</title></rect>",
                    x(visible_end.max(start)),
                    right - x(visible_end.max(start))
                );
            }
            for (reverse, scores, color) in [
                (
                    false,
                    &trace.forward,
                    if lane.kind == TssLaneKind::LocalScoreTrace {
                        "#008264"
                    } else {
                        "#1e69b9"
                    },
                ),
                (true, &trace.reverse, "#aa4669"),
            ] {
                let mut path = String::new();
                let mut hovers = String::new();
                let mut connected = false;
                for (p, score) in scores.iter().enumerate().take(visible_end).skip(start) {
                    if let Some(raw) = score.filter(|v| v.is_finite()) {
                        let value = if trace.clip_negative {
                            raw.max(0.0)
                        } else {
                            raw
                        };
                        let _ = write!(
                            path,
                            "{} {:.2},{:.2} ",
                            if connected { "L" } else { "M" },
                            x(p),
                            y(value)
                        );
                        if !connected {
                            path.push_str("l 0,0 ");
                        }
                        connected = true;
                        let details = format!(
                            "{}; {}; local strand {}; {} {}; displayed {}; footprint {}..{}; {}",
                            lane.label,
                            view.coordinate_label(p),
                            if reverse { "-" } else { "+" },
                            if lane.kind == TssLaneKind::LocalScoreTrace {
                                "computed score"
                            } else {
                                "raw"
                            },
                            raw,
                            value,
                            p + 1,
                            p + trace.motif_length_bp,
                            lane.units
                        );
                        let _ = write!(
                            hovers,
                            "<circle data-role=\"score-window\" cx=\"{:.2}\" cy=\"{:.2}\" r=\"3\" fill=\"transparent\"><title>{}</title></circle>",
                            x(p),
                            y(value),
                            xml(&details)
                        );
                    } else {
                        connected = false;
                        let _ = write!(
                            svg,
                            "<rect data-role=\"unavailable-score\" x=\"{:.2}\" y=\"{top}\" width=\"{:.2}\" height=\"120\" fill=\"#f5c878\" fill-opacity=\"0.3\"><title>Unavailable window, not zero</title></rect>",
                            x(p),
                            (x(p + 1) - x(p)).max(0.5)
                        );
                    }
                }
                let _ = write!(
                    svg,
                    "<path data-role=\"score-curve\" data-reverse=\"{reverse}\" d=\"{path}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"1.2\" stroke-linecap=\"round\" {}/>",
                    if reverse {
                        "stroke-dasharray=\"5 3\""
                    } else {
                        ""
                    }
                );
                svg.push_str(&hovers);
            }
        } else {
            let motif = matches!(lane.kind, TssLaneKind::Motif | TssLaneKind::ImportedMotif);
            let middle = motif
                || lane
                    .features
                    .iter()
                    .any(|f| f.score.is_some_and(|s| s < 0.0));
            let baseline = if middle { top + 60.0 } else { low };
            let amplitude = if middle { 60.0 } else { 120.0 };
            let color = match lane.kind {
                TssLaneKind::Structure => "#3782aa",
                TssLaneKind::Signal => "#a54164",
                TssLaneKind::Motif => "#008c73",
                _ => "#b46419",
            };
            let mut count = 0;
            for f in &lane.features {
                if f.end <= start
                    || f.start >= end
                    || (lane.kind == TssLaneKind::Motif && !f.score.is_some_and(|v| v >= 0.0))
                {
                    continue;
                }
                if f.start >= f.end || f.end > length || f.score.is_some_and(|s| !s.is_finite()) {
                    return Err("Invalid TSS feature geometry or score".into());
                }
                count += 1;
                let a = x(f.start);
                let b = x(f.end);
                let _ = write!(
                    svg,
                    "<g data-role=\"feature\" data-reverse=\"{}\"><title>{}</title>",
                    f.reverse,
                    xml(&format!(
                        "{}; {}; {}; {}",
                        f.label,
                        view.coordinate_label(f.start),
                        view.coordinate_label(f.end - 1),
                        f.details
                    ))
                );
                if motif {
                    let imported = lane.kind == TssLaneKind::ImportedMotif;
                    let value = f.score.unwrap_or(lane.scale_min);
                    let height = (if imported {
                        (value - lane.scale_min) / (lane.scale_max - lane.scale_min)
                    } else {
                        value / lane.scale_max
                    } * amplitude)
                        .max(2.0);
                    let tip = baseline + if f.reverse { height } else { -height };
                    let center = if imported { (a + b) / 2.0 } else { a };
                    let (base_a, base_b, base_y) = if imported {
                        (a, b, baseline)
                    } else {
                        (a - 3.0, a + 3.0, tip + if f.reverse { -4.0 } else { 4.0 })
                    };
                    let _ = write!(
                        svg,
                        "<path d=\"M {center:.2},{baseline:.2} V {tip:.2}\" stroke=\"{color}\"/><polygon points=\"{center:.2},{tip:.2} {base_a:.2},{base_y:.2} {base_b:.2},{base_y:.2}\" fill=\"{color}\"/>"
                    );
                } else if lane.kind == TssLaneKind::Signal {
                    if let Some(score) = f.score {
                        let value_y = baseline - score / lane.scale_max * amplitude;
                        let _ = write!(
                            svg,
                            "<rect x=\"{a:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"{color}\"/>",
                            value_y.min(baseline),
                            b - a,
                            (value_y - baseline).abs().max(1.0)
                        );
                    } else {
                        let _ = write!(
                            svg,
                            "<path data-role=\"unquantified-signal\" d=\"M {a:.2},{baseline} H {b:.2}\" stroke=\"#b47d1e\" stroke-width=\"2\"/>"
                        );
                    }
                } else {
                    let h = if f.label.contains("CDS segment") {
                        20.0
                    } else {
                        10.0
                    };
                    let _ = write!(
                        svg,
                        "<rect x=\"{a:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{h}\" fill=\"{color}\"/>",
                        top + 60.0 - h / 2.0,
                        (b - a).max(2.0)
                    );
                }
                svg.push_str("</g>");
            }
            if count == 0 {
                text(
                    &mut svg,
                    left + 8.0,
                    top + 60.0,
                    "No displayed intervals; not a measured zero",
                );
            }
        }
        svg.push_str("</g></g>");
    }
    for (n, line) in note_lines.iter().enumerate() {
        text(&mut svg, 20.0, bottom + 20.0 + n as f64 * 15.0, line);
    }
    svg.push_str("</g></svg>\n");
    if svg.len() > 32 * 1024 * 1024 {
        return Err("TSS SVG exceeds 32 MiB; select a smaller span or fewer lanes".into());
    }
    Ok(svg)
}

/// Atomically publish only a complete render. Return the exact SVG byte digest.
pub fn write_tss_view_svg(
    view: &TssSequenceView,
    options: &TssViewSvgOptions,
    path: &std::path::Path,
) -> Result<String, String> {
    use std::io::Write;
    let svg = render_tss_view_svg(view, options)?;
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(std::path::Path::new("."));
    let mut temporary = tempfile::NamedTempFile::new_in(parent).map_err(|e| e.to_string())?;
    temporary
        .write_all(svg.as_bytes())
        .map_err(|e| e.to_string())?;
    temporary.as_file().sync_all().map_err(|e| e.to_string())?;
    temporary.persist(path).map_err(|e| e.to_string())?;
    Ok(sha256_hex_bytes(svg.as_bytes()))
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reuses hand-crafted annotated sequence/profile pairs, never experimental data.
    fn fixture(minus: bool) -> (TssSequenceView, TssViewSvgOptions) {
        let (dna, report) = super::super::profile_fixture(minus);
        let view = TssSequenceView::from_dna(&dna)
            .unwrap()
            .with_profile(&report)
            .unwrap();
        let options = TssViewSvgOptions {
            start_0based: 0,
            end_0based_exclusive: 5,
            lane_indices: (0..view.lanes.len()).collect(),
            width_px: 1360,
            print_size_mm: None,
        };
        (view, options)
    }

    #[test]
    fn tss_svg_preserves_binding_scales_gaps_and_coordinates_on_both_strands() {
        for minus in [false, true] {
            let (mut view, options) = fixture(minus);
            view.title = "A&B <window>".into();
            let lane = view.lanes.iter_mut().find(|l| l.trace.is_some()).unwrap();
            lane.trace.as_mut().unwrap().forward[1] = None;
            lane.trace.as_mut().unwrap().forward[0] = Some(0.0);
            let svg = render_tss_view_svg(&view, &options).unwrap();
            assert_eq!(svg, render_tss_view_svg(&view, &options).unwrap());
            assert!(svg.contains("A&amp;B &lt;window&gt;"));
            assert!(svg.contains(&view.sequence_sha256));
            assert!(svg.contains(&sha256_hex_bytes(&serde_json::to_vec(&view).unwrap())));
            assert!(svg.contains("data-role=\"unavailable-score\""));
            assert!(svg.contains("data-role=\"terminal-unavailable\""));
            assert!(svg.contains("data-reverse=\"true\""));
            assert!(svg.contains("raw 0; displayed 0; footprint"));
            assert!(svg.contains("query complete false; truncated true"));
            assert!(svg.contains("synthetic_package_log2"));
            for pos in [0, 2, 4] {
                assert!(svg.contains(&format!(
                    ">{}</text>",
                    view.geometry.genomic_at(pos).unwrap()
                )));
            }
            let curves: Vec<_> = svg
                .split("<path data-role=\"score-curve\"")
                .skip(1)
                .collect();
            let first = curves[0].split("/>").next().unwrap();
            assert!(
                !first.contains("L "),
                "None at 1 must break the two valid windows"
            );
            assert!(first.contains("M "), "real zero must still be drawn");
            let raster = crate::svg_png::render_svg_to_png_bytes(&svg, Default::default()).unwrap();
            assert_eq!(raster.width, options.width_px);
            assert!(raster.height > 130);
            assert!(raster.bytes.starts_with(b"\x89PNG\r\n\x1a\n"));
        }
    }

    #[test]
    fn tss_svg_selection_and_atomic_export_do_not_fallback_or_overwrite_on_error() {
        let (view, mut options) = fixture(false);
        options.start_0based = 1;
        options.end_0based_exclusive = 4;
        options.lane_indices = vec![
            view.lanes
                .iter()
                .position(|l| l.kind == TssLaneKind::ImportedMotif)
                .unwrap(),
        ];
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("TSS selected view.svg");
        let hash = write_tss_view_svg(&view, &options, &path).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(hash, sha256_hex_bytes(&bytes));
        let svg = String::from_utf8(bytes.clone()).unwrap();
        assert_eq!(svg.matches("<g data-lane-id=").count(), 1);
        assert!(!svg.contains("data-role=\"score-curve\""));
        assert!(svg.contains("local 2..4"));
        options.end_0based_exclusive = 999;
        assert!(write_tss_view_svg(&view, &options, &path).is_err());
        assert_eq!(std::fs::read(&path).unwrap(), bytes);
        options.end_0based_exclusive = 4;
        options.lane_indices.push(options.lane_indices[0]);
        assert!(render_tss_view_svg(&view, &options).is_err());
        options.lane_indices.clear();
        assert!(
            !render_tss_view_svg(&view, &options)
                .unwrap()
                .contains("<g data-lane-id=")
        );
    }

    #[test]
    fn tss_svg_refuses_nonfinite_scores_and_over_budget_exports() {
        let (mut view, options) = fixture(false);
        let index = view.lanes.iter().position(|l| l.trace.is_some()).unwrap();
        view.lanes[index].trace.as_mut().unwrap().forward[0] = Some(f64::NAN);
        assert!(render_tss_view_svg(&view, &options).is_err());
        let (mut view, options) = fixture(false);
        let index = view
            .lanes
            .iter()
            .position(|l| !l.features.is_empty())
            .unwrap();
        view.lanes[index].features = vec![view.lanes[index].features[0].clone(); 10_001];
        assert!(
            render_tss_view_svg(&view, &options)
                .unwrap_err()
                .contains("paint budget")
        );
    }
}
