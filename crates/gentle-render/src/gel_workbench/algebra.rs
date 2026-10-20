//! Standalone signed density plots and complete, bounded provenance appendices.

use super::*;
use gentle_protocol::gel_workbench::{GelAlgebraPresentation, GelAlgebraVerification};

#[cfg(test)]
mod tests;

/// Project only engine-prepared density values. No CDF, QQ, significance test,
/// normalization, coefficient fitting or unrelated reference overlay is defined
/// here. Operand contributions are already rebinned/weighted by the engine.
pub fn render_algebra_svg(report: &GelAlgebraPresentation) -> String {
    let terms = report.request.as_ref().map_or(0, |r| r.terms.len().min(16));
    let mut svg = Svg::new(850 + terms * 24, "GENtle | Saved signed gel algebra");
    if report.profile_id.len() > 240 || report.label.len() > 240 {
        svg.text(
            32.0,
            36.0,
            "Unavailable: oversized algebra identity",
            18,
            80,
        );
        return svg.finish();
    }
    svg.text(32.0, 36.0, "GENtle | Saved signed gel algebra", 24, 70);
    svg.text(
        32.0,
        65.0,
        &format!("Derived result: {} | {}", report.profile_id, report.label),
        14,
        115,
    );
    svg.text(
        32.0,
        90.0,
        &format!(
            "Verification: {:?} (internal consistency only)",
            report.verification
        ),
        13,
        125,
    );
    if report.verification != GelAlgebraVerification::Verified {
        svg.text(
            32.0,
            130.0,
            report
                .reason
                .as_deref()
                .unwrap_or("No verified series; report export refused"),
            14,
            110,
        );
        return svg.finish();
    }
    let Some(series) = &report.series else {
        svg.text(
            32.0,
            130.0,
            "Unavailable: no verified density series",
            14,
            100,
        );
        return svg.finish();
    };
    if series.density.is_empty()
        || series.density.len() > MAX_BINS
        || series.weights.len() != series.density.len()
        || series.edges.len() != series.density.len() + 1
        || series.edges.iter().any(|v| !v.is_finite() || *v <= 0.0)
        || series.edges.windows(2).any(|p| p[0] >= p[1])
        || series.density.iter().flatten().any(|v| !v.is_finite())
        || series
            .weights
            .iter()
            .zip(&series.density)
            .any(|(w, d)| w.is_some() != d.is_some())
    {
        svg.text(
            32.0,
            130.0,
            "Unavailable: malformed series or 20000-bin rendering budget exceeded",
            14,
            120,
        );
        return svg.finish();
    }
    if !report.contributions.is_empty()
        && (report.contributions.len() > 16
            || report.request.as_ref().is_none_or(|request| {
                request.terms.len() > 16
                    || request
                        .terms
                        .iter()
                        .filter(|t| t.coefficient != 0.0)
                        .count()
                        != report.contributions.len()
                    || report
                        .contributions
                        .iter()
                        .zip(
                            request
                                .terms
                                .iter()
                                .enumerate()
                                .filter(|(_, t)| t.coefficient != 0.0),
                        )
                        .any(|(c, (index, t))| {
                            c.term_index != index
                                || c.profile_id != t.profile_id
                                || c.profile_id.len() > 240
                                || !c.coefficient.is_finite()
                                || c.coefficient != t.coefficient
                        })
            })
            || report.contributions.iter().any(|c| {
                c.weights.len() != series.weights.len()
                    || c.density.len() != series.density.len()
                    || c.weights
                        .iter()
                        .chain(&c.density)
                        .flatten()
                        .any(|v| !v.is_finite())
                    || c.weights
                        .iter()
                        .zip(&c.density)
                        .any(|(w, d)| w.is_some() != d.is_some())
            }))
    {
        svg.text(
            32.0,
            130.0,
            "Unavailable: malformed operand contributions",
            14,
            110,
        );
        return svg.finish();
    }
    let unit = series.size_kind.unit();
    svg.text(
        32.0,
        116.0,
        &format!(
            "Algebra input digest: {}",
            report.input_digest.as_deref().unwrap_or("unavailable")
        ),
        12,
        140,
    );
    svg.text(
        32.0,
        142.0,
        &format!(
            "Y: {} per {unit} | recorded measure: {}",
            series.signal_units, series.measure
        ),
        13,
        120,
    );
    svg.text(32.0, 166.0, "Height = saved bin weight / bin width; area retains the signed integral. No normalization or fitted scale.", 12, 135);
    svg.text(
        32.0,
        190.0,
        "Thick teal: residual. Thin styled lines: weighted operands. Black: zero; one shared axis.",
        12,
        130,
    );
    svg.text(32.0, 212.0, "Hatched: residual unavailable. Operands retain their own gaps; no line fills the residual.", 12, 130);
    let mut range = [0.0_f64, 0.0_f64];
    for value in series
        .density
        .iter()
        .chain(report.contributions.iter().flat_map(|c| c.density.iter()))
        .flatten()
    {
        range[0] = range[0].min(*value);
        range[1] = range[1].max(*value);
    }
    if !(range[1] - range[0]).is_finite() {
        svg.text(
            32.0,
            232.0,
            "Unavailable: finite density range cannot be rendered",
            14,
            100,
        );
        return svg.finish();
    }
    let fallback_axis = range[0] == range[1];
    if fallback_axis {
        range = [-1.0, 1.0];
    }
    let plot = Plot {
        x: 105.0,
        y: 270.0,
        w: 908.0,
        h: 300.0,
        lo: series.edges[0],
        hi: *series.edges.last().unwrap(),
    };
    plot.axes(&mut svg, range, &format!("Size ({unit})"));
    for (i, value) in series.density.iter().enumerate() {
        if value.is_none() {
            svg.rect(
                plot.x(series.edges[i]),
                plot.y,
                plot.x(series.edges[i + 1]) - plot.x(series.edges[i]),
                plot.h,
                "url(#missing)",
            );
        }
    }
    let zero = plot.y(0.0, range);
    let _ = writeln!(
        svg.0,
        "<path id=\"algebra-zero\" d=\"M{} {zero}L{} {zero}\" stroke=\"#23313d\" fill=\"none\"><title>Zero signed signal</title></path>",
        plot.x,
        plot.x + plot.w
    );
    svg.text(plot.x + plot.w + 8.0, zero + 4.0, "0 (zero)", 12, 12);
    for c in &report.contributions {
        let (color, dash) = contribution_style(c.term_index);
        let path = step_path(plot, range, &series.edges, &c.density);
        let _ = writeln!(
            svg.0,
            "<path id=\"algebra-operand-{}\" d=\"{path}\" stroke=\"{color}\" stroke-width=\"1.2\" stroke-dasharray=\"{dash}\" fill=\"none\"><title>Term {}: {} x {}; coefficient-weighted integral per {unit}; gaps are unavailable</title></path>",
            c.term_index,
            c.term_index + 1,
            c.coefficient,
            xml(&c.profile_id)
        );
    }
    // Paint the residual last so cancellation to zero cannot be hidden by terms.
    let path = step_path(plot, range, &series.edges, &series.density);
    svg.path("algebra-density", &path, LEFT, false);
    if series.density.iter().all(Option::is_none) {
        svg.text(
            105.0,
            242.0,
            if fallback_axis {
                "No evaluated bins; +/-1 is a display fallback, NOT measured zero signal"
            } else {
                "No evaluated bins in residual; visible operands do not establish residual coverage"
            },
            12,
            120,
        );
    } else if series.density.iter().flatten().all(|v| *v == 0.0) {
        svg.text(
            105.0,
            242.0,
            if fallback_axis {
                "All evaluated bins are exactly zero; display range +/-1 (gaps remain unavailable)"
            } else {
                "All evaluated bins are exactly zero; shared axis includes weighted operands (gaps stay unavailable)"
            },
            12,
            120,
        );
    }
    svg.text(
        32.0,
        615.0,
        "Weighted operand legend (same axis; no fitting; full identities/calibrations in provenance pages)",
        13,
        125,
    );
    if let Some(request) = &report.request {
        for (i, term) in request.terms.iter().take(16).enumerate() {
            let contribution = report.contributions.iter().find(|c| c.term_index == i);
            if contribution.is_some() {
                let (color, dash) = contribution_style(i);
                let y = 641.0 + i as f64 * 24.0 - 4.0;
                let _ = writeln!(
                    svg.0,
                    "<path id=\"algebra-operand-legend-{i}\" d=\"M32 {y}L65 {y}\" stroke=\"{color}\" stroke-width=\"1.2\" stroke-dasharray=\"{dash}\"/>"
                );
            }
            let identity = report
                .inputs
                .iter()
                .find(|p| p.id == term.profile_id)
                .map(|p| {
                    if p.evidence_kind == GelEvidenceKind::SignedResidual {
                        "derived operand; see dependency receipt".into()
                    } else {
                        format!("run {} / lane {}", p.run_id, p.lane_id)
                    }
                })
                .unwrap_or_else(|| "missing operand identity".into());
            svg.text(
                76.0,
                641.0 + i as f64 * 24.0,
                &format!(
                    "Term {}: {} x {} | {identity}{}",
                    i + 1,
                    term.coefficient,
                    term.profile_id,
                    if term.coefficient == 0.0 {
                        " | zero coefficient: provenance/grid only; not drawn"
                    } else if contribution.is_none() {
                        " | overlay unavailable; re-inspect legacy presentation"
                    } else {
                        ""
                    }
                ),
                12,
                140,
            );
        }
    }
    let y = 680.0 + terms as f64 * 24.0;
    svg.text(32.0, y, "Exploratory signed arithmetic, not gene-specific enrichment or a validated additive reaction model.", 12, 135);
    svg.text(32.0, y + 25.0, "Image-derived inputs remain display diagnostics, not fluorescence, concentration, molecule counts or yield.", 12, 135);
    svg.text(
        32.0,
        y + 50.0,
        "Every operand's processing, warnings and hashes are retained in the appendix and JSON.",
        12,
        135,
    );
    svg.finish()
}

fn contribution_style(index: usize) -> (&'static str, &'static str) {
    let colors = ["#3156a3", "#bd5b21", "#8b518b", "#5f7424"];
    let dashes = ["none", "7 3", "2 3", "8 3 2 3"];
    (colors[index % 4], dashes[(index / 4) % 4])
}

fn step_path(plot: Plot, range: [f64; 2], edges: &[f64], values: &[Option<f64>]) -> String {
    let mut path = String::new();
    let mut connected = false;
    for (i, value) in values.iter().enumerate() {
        if let Some(value) = value {
            let _ = write!(
                path,
                "{} {:.5} {:.5} L {:.5} {:.5} ",
                if connected { "L" } else { "M" },
                plot.x(edges[i]),
                plot.y(*value, range),
                plot.x(edges[i + 1]),
                plot.y(*value, range)
            );
            connected = true;
        } else {
            connected = false;
        }
    }
    path
}

/// Append complete wrapped metadata, or refuse explicitly rather than silently
/// omitting provenance. The JSON export is the fallback for oversized metadata.
pub fn render_algebra_provenance_svgs(
    report: &GelAlgebraPresentation,
) -> Result<Vec<String>, String> {
    if report.verification != GelAlgebraVerification::Verified || report.inputs.len() > 64 {
        return Err("Algebra provenance requires a verified, bounded presentation".into());
    }
    if report.inputs.iter().any(|p| {
        p.id.len() > 240
            || p.label.len() > 240
            || p.run_id.len() > 240
            || p.lane_id.len() > 240
            || p.source_sha256.len() != 64
            || p.signal_units.len() > 240
            || p.measure.len() > 240
            || p.weights.len() > MAX_BINS
            || p.processing.len() > 128
            || p.warnings.len() > 128
            || p.processing
                .iter()
                .chain(&p.warnings)
                .any(|s| s.len() > 2048)
    }) || report.warnings.len() > 8328
        || report.warnings.iter().any(|s| s.len() > 2048)
    {
        return Err("Algebra provenance metadata exceeds rendering bounds".into());
    }
    let mut lines = Vec::new();
    let mut add = |value: String| -> Result<(), String> {
        let chars: Vec<_> = value.chars().collect();
        let mut start = 0;
        while start < chars.len() {
            let mut end = (start + 78).min(chars.len());
            if end < chars.len()
                && let Some(space) = chars[start..end].iter().rposition(|c| c.is_whitespace())
                && space > 0
            {
                end = start + space + 1;
            }
            lines.push(chars[start..end].iter().collect::<String>());
            start = end;
            if lines.len() > 50 * 64 {
                return Err(
                    "Algebra provenance exceeds 64 pages; export JSON or reduce selected results"
                        .into(),
                );
            }
        }
        Ok(())
    };
    add(format!("Derived result: {}", report.profile_id))?;
    add("Algebra input digest:".into())?;
    add(report
        .input_digest
        .as_deref()
        .unwrap_or("unavailable")
        .into())?;
    if let Some(request) = &report.request {
        for term in &request.terms {
            add(format!(
                "Coefficient {} x {}{}",
                term.coefficient,
                term.profile_id,
                if term.coefficient == 0.0 {
                    " (grid/provenance only; no signal or gap exclusion)"
                } else {
                    ""
                }
            ))?;
        }
    }
    for profile in &report.inputs {
        add(format!(
            "Operand/dependency: {} | {}",
            profile.id, profile.label
        ))?;
        if profile.evidence_kind == GelEvidenceKind::SignedResidual {
            add("Derived dependency: inherited run/lane fields are not its identity".into())?;
            add("Algebra input digest:".into())?;
            add(profile.source_sha256.clone())?;
        } else {
            add(format!(
                "Run {} | lane {} | evidence {:?}",
                profile.run_id, profile.lane_id, profile.evidence_kind
            ))?;
            add("Source SHA-256:".into())?;
            add(profile.source_sha256.clone())?;
        }
        add(format!(
            "Units {} | measure {} | {} bins, {} unavailable",
            profile.signal_units,
            profile.measure,
            profile.weights.len(),
            profile.weights.iter().filter(|w| w.is_none()).count()
        ))?;
        for value in &profile.processing {
            add(format!("Processing/calibration: {value}"))?;
        }
        for value in &profile.warnings {
            add(format!("Operand warning: {value}"))?;
        }
    }
    for warning in &report.warnings {
        add(format!("Warning: {warning}"))?;
    }
    Ok(lines
        .chunks(50)
        .enumerate()
        .map(|(i, lines)| {
            let mut svg = Svg::new(1020, "GENtle | Algebra provenance");
            svg.text(
                32.0,
                36.0,
                &format!("Algebra provenance | page {}", i + 1),
                22,
                75,
            );
            for (j, line) in lines.iter().enumerate() {
                svg.text(32.0, 76.0 + j as f64 * 18.0, line, 12, 140);
            }
            svg.finish()
        })
        .collect())
}
