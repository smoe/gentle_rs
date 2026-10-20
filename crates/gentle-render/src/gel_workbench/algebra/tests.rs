//! Synthetic frozen presentations, not detector measurements or engine replay.

use super::*;
use gentle_protocol::{gel_image::GelSizeKind, gel_workbench::*};

fn presentation() -> GelAlgebraPresentation {
    GelAlgebraPresentation {
        schema: "gentle.gel_algebra_presentation.v1".into(),
        profile_id: "synthetic residual".into(),
        label: "Saved result".into(),
        verification: GelAlgebraVerification::Verified,
        reason: None,
        input_digest: Some("d".repeat(64)),
        request: Some(GelAlgebraRequest {
            output_id: "synthetic residual".into(),
            range: [100.0, 300.0],
            terms: vec![GelAlgebraTerm {
                profile_id: "a".into(),
                coefficient: 1.0,
            }],
        }),
        series: Some(GelAlgebraSeries {
            size_kind: GelSizeKind::LinearDnaBp,
            signal_units: "synthetic_area".into(),
            measure: "synthetic_integral".into(),
            edges: vec![100.0, 150.0, 200.0, 300.0],
            weights: vec![Some(-1.0), None, Some(3.0)],
            density: vec![Some(-0.02), None, Some(0.03)],
        }),
        contributions: vec![],
        reference: None,
        inputs: vec![GelProfile {
            id: "a".into(),
            label: "Synthetic input".into(),
            run_id: "source run".into(),
            lane_id: "C1".into(),
            source_sha256: "a".repeat(64),
            evidence_kind: GelEvidenceKind::ImageDisplay,
            size_kind: GelSizeKind::LinearDnaBp,
            signal_units: "synthetic_area".into(),
            measure: "synthetic_integral".into(),
            edges: vec![100.0, 300.0],
            weights: vec![Some(2.0)],
            processing: vec!["Calibration L2; background 0; masked ranges none".into()],
            warnings: vec!["Display-only source; not calibrated detector data".into()],
        }],
        warnings: vec!["No biological inference".into()],
    }
}

fn reference_presentation() -> GelAlgebraPresentation {
    let mut r = presentation();
    let mut p = r.inputs[0].clone();
    p.id = "custom<&".into();
    p.signal_units = "different_reference_area".into();
    p.edges = vec![90.0, 125.0, 175.0, 250.0, 350.0];
    p.weights = vec![Some(70.0), None, Some(150.0), Some(0.0)];
    r.reference = Some(GelAlgebraReferencePanel {
        series: GelAlgebraSeries {
            size_kind: p.size_kind,
            signal_units: p.signal_units.clone(),
            measure: p.measure.clone(),
            edges: p.edges.clone(),
            weights: p.weights.clone(),
            density: vec![Some(2.0), None, Some(2.0), Some(0.0)],
        },
        profile: p,
    });
    r
}

#[test]
fn signed_algebra_reference_has_native_bins_an_independent_axis_and_no_residual_rescale() {
    let plain = render_algebra_svg(&presentation());
    let with = reference_presentation();
    let svg = render_algebra_svg(&with);
    let residual_path = |text: &str| {
        text.split("id=\"algebra-density\" d=\"")
            .nth(1)
            .unwrap()
            .split('"')
            .next()
            .unwrap()
            .to_string()
    };
    assert_eq!(residual_path(&plain), residual_path(&svg));
    let path = svg
        .split("id=\"algebra-reference-density\" d=\"")
        .nth(1)
        .unwrap()
        .split('"')
        .next()
        .unwrap();
    assert_eq!(
        path.matches('M').count(),
        2,
        "native missing bin breaks the reference path"
    );
    assert!(
        path.starts_with("M 105.00000"),
        "native 90..125 bin is only display-clipped at 100"
    );
    assert!(svg.contains("UNSCALED; separate Y axis"));
    assert!(svg.contains("Reference: custom&lt;&amp;"));
    assert!(svg.contains("different_reference_area per bp"));
    assert!(svg.contains("Compare positions, not amplitudes"));
    assert!(svg.contains("registration remain undecided"));
    assert!(!svg.contains("id=\"cdf") && !svg.contains("id=\"qq"));
    let appendix = render_algebra_provenance_svgs(&with).unwrap().join("\n");
    assert!(appendix.contains("NOT an operand"));
    assert!(appendix.contains("Calibration L2"));
    assert!(appendix.contains(&with.reference.as_ref().unwrap().profile.source_sha256));
    assert!(!plain.contains("algebra-reference-panel"));
}

#[test]
fn signed_algebra_reference_preserves_zero_missing_and_refuses_malformed_series() {
    for missing in [false, true] {
        let mut r = reference_presentation();
        let reference = r.reference.as_mut().unwrap();
        let value = if missing { None } else { Some(0.0) };
        reference.profile.weights.fill(value);
        reference.series.weights.fill(value);
        reference.series.density.fill(value);
        let svg = render_algebra_svg(&r);
        assert!(svg.contains(if missing {
            "No evaluated reference bins"
        } else {
            "All evaluated reference bins are zero"
        }));
        assert!(svg.contains("display fallback"));
    }
    for mode in 0..3 {
        let mut r = reference_presentation();
        let reference = r.reference.as_mut().unwrap();
        match mode {
            0 => reference.series.density[0] = Some(f64::NAN),
            1 => reference.series.weights[0] = Some(999.0),
            _ => reference.profile.evidence_kind = GelEvidenceKind::SignedResidual,
        }
        let svg = render_algebra_svg(&r);
        assert!(svg.contains("malformed custom-only reference"));
        assert!(!svg.contains("id=\"algebra-density\""));
        assert!(!svg.contains("id=\"algebra-reference-density\""));
    }
}

#[test]
fn signed_algebra_weighted_operands_share_axis_keep_gaps_and_use_ordered_styles() {
    let mut r = presentation();
    r.request.as_mut().unwrap().terms = vec![
        GelAlgebraTerm {
            profile_id: "a<&".into(),
            coefficient: 2.0,
        },
        GelAlgebraTerm {
            profile_id: "b".into(),
            coefficient: -0.5,
        },
        GelAlgebraTerm {
            profile_id: "grid-only".into(),
            coefficient: 0.0,
        },
    ];
    r.contributions = vec![
        GelAlgebraContribution {
            term_index: 0,
            profile_id: "a<&".into(),
            coefficient: 2.0,
            weights: vec![Some(10.0), None, Some(20.0)],
            density: vec![Some(0.2), None, Some(0.2)],
        },
        GelAlgebraContribution {
            term_index: 1,
            profile_id: "b".into(),
            coefficient: -0.5,
            weights: vec![Some(-15.0), Some(-15.0), Some(-30.0)],
            density: vec![Some(-0.3), Some(-0.3), Some(-0.3)],
        },
    ];
    r.series.as_mut().unwrap().weights = vec![Some(-5.0), None, Some(-10.0)];
    r.series.as_mut().unwrap().density = vec![Some(-0.1), None, Some(-0.1)];
    let svg = render_algebra_svg(&r);
    for (id, moves) in [("algebra-operand-0", 2), ("algebra-operand-1", 1)] {
        let path = svg
            .split(&format!("id=\"{id}\" d=\""))
            .nth(1)
            .unwrap()
            .split('"')
            .next()
            .unwrap();
        assert_eq!(path.matches('M').count(), moves);
    }
    assert!(svg.contains("Term 1: 2 x a&lt;&amp;"));
    assert!(svg.contains("Term 2: -0.5 x b"));
    assert!(svg.contains("zero coefficient: provenance/grid only; not drawn"));
    assert!(!svg.contains("id=\"algebra-operand-2\""));
    assert!(svg.contains("0.2"));
    assert!(svg.contains("-0.3"));
    // Shared range [-0.3, 0.2], plot y=270..570: positive/negative terms
    // reach opposite ends; the residual is inside it, not independently fitted.
    assert!(svg.contains("M 105.00000 270.00000"));
    assert!(svg.contains("M 105.00000 570.00000"));
    assert!(
        svg.find("id=\"algebra-density\"").unwrap() > svg.find("id=\"algebra-operand-1\"").unwrap()
    );
    assert!(!svg.contains("id=\"cdf"));
    assert!(!svg.contains("id=\"qq"));
    assert_eq!(svg, render_algebra_svg(&r));
    for index in 0..16 {
        for previous in 0..index {
            assert_ne!(contribution_style(index), contribution_style(previous));
        }
    }
}

#[test]
fn signed_algebra_malformed_operand_projection_never_draws_a_curve() {
    for fault in 0..6 {
        let mut r = presentation();
        r.contributions = vec![GelAlgebraContribution {
            term_index: 0,
            profile_id: "a".into(),
            coefficient: 1.0,
            weights: r.series.as_ref().unwrap().weights.clone(),
            density: r.series.as_ref().unwrap().density.clone(),
        }];
        match fault {
            0 => r.contributions[0].density.pop().map(|_| ()).unwrap(),
            1 => r.contributions[0].coefficient = -1.0,
            2 => r.contributions[0].term_index = 1,
            3 => r.contributions[0].density[0] = Some(f64::INFINITY),
            4 => r.contributions[0].weights[0] = None,
            _ => r.contributions.push(r.contributions[0].clone()),
        }
        let svg = render_algebra_svg(&r);
        assert!(svg.contains("malformed operand contributions"));
        assert!(!svg.contains("id=\"algebra-density\""));
        assert!(!svg.contains("id=\"algebra-operand-0\""));
    }
}

#[test]
fn signed_algebra_plot_has_zero_density_units_gaps_and_no_distribution_panels() {
    let r = presentation();
    let svg = render_algebra_svg(&r);
    for text in [
        "algebra-zero",
        "0 (zero)",
        "algebra-density",
        "per bp",
        "url(#missing)",
        "Algebra input digest",
        "No normalization",
        "not gene-specific",
        "-0.02",
    ] {
        assert!(svg.contains(text), "{text}");
    }
    assert!(!svg.contains("id=\"cdf"));
    assert!(!svg.contains("id=\"qq"));
    assert_eq!(svg, render_algebra_svg(&r));
}

#[test]
fn signed_algebra_gaps_split_paths_and_all_missing_is_not_zero() {
    let mut r = presentation();
    let svg = render_algebra_svg(&r);
    let path = svg
        .split("id=\"algebra-density\" d=\"")
        .nth(1)
        .unwrap()
        .split('"')
        .next()
        .unwrap();
    assert_eq!(path.matches('M').count(), 2);
    r.series.as_mut().unwrap().density = vec![None; 3];
    r.series.as_mut().unwrap().weights = vec![None; 3];
    assert!(render_algebra_svg(&r).contains("No evaluated bins"));
    assert!(render_algebra_svg(&r).contains("display fallback, NOT measured zero"));
    r.series.as_mut().unwrap().density = vec![Some(0.0); 3];
    r.series.as_mut().unwrap().weights = vec![Some(0.0); 3];
    assert!(render_algebra_svg(&r).contains("All evaluated bins are exactly zero"));
    r.series.as_mut().unwrap().size_kind = GelSizeKind::SdsProteinKda;
    assert!(render_algebra_svg(&r).contains("per kDa"));
}

#[test]
fn signed_algebra_provenance_retains_all_calibrations_hashes_and_escaped_labels() {
    let mut r = presentation();
    let mut second = r.inputs[0].clone();
    second.id = "b<script>".into();
    second.source_sha256 = "b".repeat(64);
    second.processing = vec!["Calibration L1 & reviewed <background>".into()];
    second.warnings = vec!["Second-operand warning".into()];
    r.inputs.push(second);
    let warning = "Every operand retains its own calibration/processing. Different ladders or registration uncertainty require review; subtraction does not align them.";
    r.warnings.push(warning.into());
    let pages = render_algebra_provenance_svgs(&r).unwrap().join("\n");
    for text in [
        "a".repeat(64),
        "b".repeat(64),
        "Calibration L2".into(),
        "Calibration L1 &amp; reviewed &lt;background&gt;".into(),
        "Second-operand warning".into(),
    ] {
        assert!(pages.contains(&text), "{text}");
    }
    for word in warning.split_whitespace() {
        assert!(pages.contains(word), "Provenance split a word: {word}");
    }
    assert!(!pages.contains(">.</text>"));
    assert!(!pages.contains("<script>"));
    r.inputs[0].evidence_kind = GelEvidenceKind::SignedResidual;
    r.inputs[0].run_id = "not-a-derived-run".into();
    let pages = render_algebra_provenance_svgs(&r).unwrap().join("\n");
    assert!(pages.contains("inherited run/lane fields are not its identity"));
    assert!(!pages.contains("not-a-derived-run"));
}

#[test]
fn signed_algebra_unverified_and_malformed_series_never_become_a_plausible_curve() {
    let mut r = presentation();
    r.verification = GelAlgebraVerification::Mismatch;
    assert!(!render_algebra_svg(&r).contains("id=\"algebra-density\""));
    assert!(render_algebra_provenance_svgs(&r).is_err());
    r = presentation();
    r.series.as_mut().unwrap().edges[1] = f64::NAN;
    assert!(render_algebra_svg(&r).contains("malformed series"));
    r = presentation();
    r.series.as_mut().unwrap().density.resize(20_001, Some(1.0));
    assert!(render_algebra_svg(&r).contains("20000-bin"));
}

#[test]
fn signed_algebra_full_bin_budget_and_provenance_overflow_are_explicit() {
    let mut r = presentation();
    let s = r.series.as_mut().unwrap();
    s.edges = (0..=20_000).map(|i| 100.0 + i as f64).collect();
    s.weights = vec![Some(-0.01); 20_000];
    s.density = s.weights.clone();
    let svg = render_algebra_svg(&r);
    let path = svg
        .split("id=\"algebra-density\" d=\"")
        .nth(1)
        .unwrap()
        .split('"')
        .next()
        .unwrap();
    assert_eq!(path.matches('L').count(), 39_999);
    r.warnings = vec!["w".repeat(2048); 200];
    assert!(
        render_algebra_provenance_svgs(&r)
            .unwrap_err()
            .contains("64 pages")
    );
}
