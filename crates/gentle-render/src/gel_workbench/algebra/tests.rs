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
