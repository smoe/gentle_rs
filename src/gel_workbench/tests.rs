//! Independently calculable synthetic profiles and images; no detector samples
//! or TapeStation region means are manufactured from a rendered gel.

use super::*;

fn profile(id: &str, edges: &[f64], weights: &[Option<f64>]) -> GelProfile {
    GelProfile {
        id: id.into(),
        label: id.into(),
        run_id: "synthetic".into(),
        lane_id: id.into(),
        source_sha256: "a".repeat(64),
        evidence_kind: GelEvidenceKind::VendorTrace,
        size_kind: GelSizeKind::LinearDnaBp,
        signal_units: "synthetic_integrated_signal".into(),
        measure: "synthetic_bin_area".into(),
        edges: edges.into(),
        weights: weights.into(),
        processing: vec!["hand-calculated synthetic bin integrals".into()],
        warnings: vec![],
    }
}
fn request() -> GelCompareRequest {
    GelCompareRequest {
        comparison_id: "compare".into(),
        left_profile: "a".into(),
        right_profile: "b".into(),
        range: [100.0, 300.0],
        short_threshold: Some(200.0),
    }
}

#[test]
fn identical_weighted_profiles_have_identity_qq_and_no_inferential_p_value() {
    let a = profile("a", &[100.0, 200.0, 300.0], &[Some(1.0), Some(1.0)]);
    let mut b = a.clone();
    b.id = "b".into();
    let r = compare(&a, &b, &request()).unwrap();
    assert_eq!(r.ks_distance, 0.0);
    assert_eq!(r.wasserstein_distance, 0.0);
    assert_eq!(r.p_value, None);
    assert_eq!(r.left_summary.median, 200.0);
    assert_eq!(r.left_summary.fraction_below_threshold, Some(0.5));
    assert_eq!(r.qq.len(), 99);
    assert!(r.qq.iter().all(|p| p.left == p.right));
}

#[test]
fn shifted_uniform_bins_have_known_ks_wasserstein_and_qq() {
    let a = profile("a", &[100.0, 200.0, 300.0], &[Some(1.0), Some(0.0)]);
    let b = profile("b", &[100.0, 200.0, 300.0], &[Some(0.0), Some(1.0)]);
    let r = compare(&a, &b, &request()).unwrap();
    assert_eq!(r.ks_distance, 1.0);
    assert!((r.wasserstein_distance - 100.0).abs() < 1e-10);
    assert_eq!(r.left_summary.median, 150.0);
    assert_eq!(r.right_summary.median, 250.0);
    assert!(
        r.qq.iter()
            .all(|q| ((q.right - q.left) - 100.0).abs() < 1e-10)
    );
}

#[test]
fn rebin_preserves_integrals_and_excludes_missing_support_for_both_profiles() {
    let a = profile("a", &[100.0, 300.0], &[Some(4.0)]);
    let b = profile(
        "b",
        &[100.0, 150.0, 200.0, 250.0, 300.0],
        &[Some(1.0), None, Some(1.0), Some(1.0)],
    );
    let r = compare(&a, &b, &request()).unwrap();
    assert_eq!(r.left_weights, vec![Some(1.0), None, Some(1.0), Some(1.0)]);
    assert_eq!(r.excluded_intervals, vec![[150.0, 200.0]]);
    assert_eq!(r.left_summary.total_weight, 3.0);
    assert_eq!(r.ks_distance, 0.0);
    assert_eq!(r.wasserstein_distance, 0.0);
}

#[test]
fn signed_algebra_retains_negatives_gaps_and_explicit_coefficients() {
    let a = profile("a", &[100.0, 200.0, 300.0], &[Some(1.0), Some(3.0)]);
    let b = profile("b", &[100.0, 200.0, 300.0], &[Some(2.0), None]);
    let request = GelAlgebraRequest {
        output_id: "residual".into(),
        terms: vec![
            GelAlgebraTerm {
                profile_id: "a".into(),
                coefficient: 1.0,
            },
            GelAlgebraTerm {
                profile_id: "b".into(),
                coefficient: -1.0,
            },
        ],
        range: [100.0, 300.0],
    };
    let r = algebra(&[&a, &b], &request).unwrap();
    assert_eq!(r.weights, vec![Some(-1.0), None]);
    assert_eq!(r.evidence_kind, GelEvidenceKind::SignedResidual);
    let mut b = r.clone();
    b.id = "b".into();
    let mut a = r.clone();
    a.id = "a".into();
    assert!(compare(&a, &b, &super::tests::request()).is_err());
}

#[test]
fn refuses_incompatible_measures_negative_probabilities_and_empty_distributions() {
    let mut a = profile("a", &[100.0, 300.0], &[Some(1.0)]);
    let mut b = a.clone();
    b.id = "b".into();
    b.evidence_kind = GelEvidenceKind::ImageDisplay;
    assert!(compare(&a, &b, &request()).is_err());
    b = a.clone();
    b.id = "b".into();
    b.measure = "other integration measure".into();
    assert!(compare(&a, &b, &request()).is_err());
    a.weights[0] = Some(-1.0);
    assert!(validate_profile(&a).is_err());
    a.weights[0] = Some(0.0);
    b = a.clone();
    b.id = "b".into();
    assert!(compare(&a, &b, &request()).is_err());
    a.edges[1] = a.edges[0];
    assert!(validate_profile(&a).is_err());
}

#[test]
fn automatic_proposals_find_synthetic_bands_without_assigning_ladder_identities() {
    let mut pixels = image::GrayImage::from_pixel(30, 100, image::Luma([255]));
    for y in [20, 70] {
        for x in 5..25 {
            pixels.put_pixel(x, y, image::Luma([0]));
        }
    }
    let mut bytes = std::io::Cursor::new(Vec::new());
    image::DynamicImage::ImageLuma8(pixels)
        .write_to(&mut bytes, image::ImageFormat::Png)
        .unwrap();
    let image =
        crate::gel_image::import_gel_image_bytes("image", "synthetic.png", bytes.get_ref(), None)
            .unwrap();
    let request = GelDetectRequest {
        image_id: "image".into(),
        lanes: vec![GelImageLane {
            id: "sample".into(),
            label: "sample".into(),
            min: GelImagePoint { x: 5.0, y: 0.0 },
            max: GelImagePoint { x: 24.0, y: 99.0 },
        }],
        migration: GelMigrationDirection::Down,
        dark_bands: true,
        minimum_prominence: 0.5,
        minimum_distance_px: 5,
        smoothing_radius_px: 0,
    };
    let r = detect_bands(&image, &request).unwrap();
    assert_eq!(r.bands.len(), 2);
    assert_eq!(r.bands[0].center.y, 20.0);
    assert_eq!(r.bands[1].center.y, 70.0);
    assert!(r.bands.iter().all(|b| b.label.contains("unconfirmed")));
}

#[test]
fn shared_shell_persists_profiles_comparison_and_algebra_without_export_side_effects() {
    use crate::{
        engine::{Engine, GentleEngine, Operation},
        engine_shell::{execute_shell_command, parse_shell_line},
    };
    let mut engine = GentleEngine::new();
    for p in [
        profile("a", &[100.0, 200.0, 300.0], &[Some(1.0), Some(1.0)]),
        profile("b", &[100.0, 200.0, 300.0], &[Some(0.0), Some(2.0)]),
    ] {
        let request = GelWorkbenchRequest::ImportProfile { profile: p };
        let line = format!(
            "gel-image workbench {}",
            crate::engine_shell::quote_shell_arg(&serde_json::to_string(&request).unwrap())
        );
        let command = parse_shell_line(&line).unwrap();
        assert!(command.is_state_mutating());
        let r = execute_shell_command(&mut engine, &command).unwrap();
        assert!(r.state_changed);
    }
    let r = engine
        .apply(Operation::GelWorkbench {
            request: Box::new(GelWorkbenchRequest::Compare(request())),
        })
        .unwrap();
    assert!(r.gel_workbench.unwrap().comparison.is_some());
    let snapshot = serde_json::to_value(engine.state()).unwrap();
    assert_eq!(
        snapshot["gel_images"]["workbench"]["comparisons"]["compare"]["p_value"],
        serde_json::Value::Null
    );
    let before = engine.mutation_revision();
    engine
        .apply(Operation::GelWorkbench {
            request: Box::new(GelWorkbenchRequest::Inspect {}),
        })
        .unwrap();
    assert_eq!(before, engine.mutation_revision());
}

#[test]
fn numeric_csv_retains_exact_lf_crlf_bytes_gaps_and_declared_measure() {
    let request = GelProfileFileRequest {
        path: "trace.csv".into(),
        profile_id: "csv".into(),
        label: "CSV".into(),
        run_id: "run".into(),
        lane_id: "sample".into(),
        size_kind: GelSizeKind::LinearDnaBp,
        signal_units: "vendor_area".into(),
        measure: "declared_vendor_bin_integrals".into(),
    };
    for newline in ["\n", "\r\n"] {
        let csv = ["left_edge,right_edge,weight", "100,150,2", "200,300,"].join(newline) + newline;
        let (p, source) = import_profile_csv(&request, csv.as_bytes()).unwrap();
        assert_eq!(p.edges, vec![100.0, 150.0, 200.0, 300.0]);
        assert_eq!(p.weights, vec![Some(2.0), None, None]);
        assert_eq!(source.sha256, sha256_hex_bytes(csv.as_bytes()));
        assert_eq!(p.source_sha256, source.sha256);
        assert_eq!(
            STANDARD.decode(source.original_base64).unwrap(),
            csv.as_bytes()
        );
    }
    for bad in [
        "size,fluorescence\n100,8\n",
        "left_edge,right_edge,weight\n100,200,-1\n",
        "left_edge,right_edge,weight\n100,200,1\n150,250,2\n",
    ] {
        assert!(import_profile_csv(&request, bad.as_bytes()).is_err());
    }
}

#[test]
fn saved_comparison_json_export_is_read_only_no_clobber_and_tamper_checked() {
    use crate::engine::{Engine, GentleEngine, Operation};
    let mut engine = GentleEngine::new();
    for p in [
        profile("a", &[100.0, 300.0], &[Some(2.0)]),
        profile("b", &[100.0, 300.0], &[Some(3.0)]),
    ] {
        engine
            .apply(Operation::GelWorkbench {
                request: Box::new(GelWorkbenchRequest::ImportProfile { profile: p }),
            })
            .unwrap();
    }
    engine
        .apply(Operation::GelWorkbench {
            request: Box::new(GelWorkbenchRequest::Compare(request())),
        })
        .unwrap();
    let mut reopened = GentleEngine::from_state(
        serde_json::from_value(serde_json::to_value(engine.state()).unwrap()).unwrap(),
    );
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("saved comparison's report.json");
    let request = GelWorkbenchExportRequest {
        path: path.display().to_string(),
        title: "Synthetic before/after".into(),
        analysis_ids: vec![],
        comparison_ids: vec!["compare".into()],
    };
    let op = Operation::GelWorkbench {
        request: Box::new(GelWorkbenchRequest::Export(request.clone())),
    };
    let before = serde_json::to_value(reopened.state()).unwrap();
    let r = reopened.apply(op.clone()).unwrap();
    assert_eq!(before, serde_json::to_value(reopened.state()).unwrap());
    let bytes = std::fs::read(&path).unwrap();
    assert_eq!(
        r.gel_workbench.unwrap().exported_sha256,
        Some(sha256_hex_bytes(&bytes))
    );
    assert!(reopened.apply(op).is_err());
    assert_eq!(bytes, std::fs::read(&path).unwrap());
    let mut tampered = serde_json::to_value(reopened.state()).unwrap();
    tampered["gel_images"]["workbench"]["comparisons"]["compare"]["p_value"] =
        serde_json::json!(0.001);
    let mut altered = GentleEngine::from_state(serde_json::from_value(tampered).unwrap());
    let other = dir.path().join("tampered.json");
    let mut request = request;
    request.path = other.display().to_string();
    assert!(
        altered
            .apply(Operation::GelWorkbench {
                request: Box::new(GelWorkbenchRequest::Export(request))
            })
            .is_err()
    );
    assert!(!other.exists());
}

#[test]
fn operation_source_summary_never_exposes_native_archive_payloads() {
    let source = GelSourceRecord {
        request: GelSourceImportRequest {
            image_id: "native".into(),
            path: "private.D5000".into(),
            pdf_object: None,
            allow_native_thumbnail: true,
        },
        image_id: "native".into(),
        source_name: "private.D5000".into(),
        source_sha256: "a".repeat(64),
        extracted_image_sha256: "b".repeat(64),
        original_base64: STANDARD.encode(b"private operational logs"),
        extraction: "thumbnail_only".into(),
        warnings: vec![],
    };
    let json = serde_json::to_value(GelSourceSummary::from(&source)).unwrap();
    assert!(json.get("original_base64").is_none());
    assert_eq!(json["source_sha256"], source.source_sha256);
}

#[test]
fn wasserstein_large_finite_axis_does_not_overflow_intermediate_product() {
    let edges = [1.0, 2.0, 1e308, 1.5e308];
    let a = profile("a", &edges, &[Some(1.0), Some(0.0), Some(0.0)]);
    let b = profile("b", &edges, &[Some(0.0), Some(0.0), Some(1.0)]);
    let r = compare(
        &a,
        &b,
        &GelCompareRequest {
            range: [1.0, 1.5e308],
            short_threshold: None,
            ..request()
        },
    )
    .unwrap();
    assert!(r.wasserstein_distance.is_finite());
    assert!((r.wasserstein_distance / 1e308 - 1.25).abs() < 1e-12);
}

#[test]
fn native_alpha_is_composited_without_discarding_sixteen_bit_source() {
    let mut pixels = image::ImageBuffer::<image::Rgba<u16>, Vec<u16>>::from_pixel(
        2,
        2,
        image::Rgba([0, 0, 0, 0]),
    );
    pixels.put_pixel(1, 1, image::Rgba([0, 0, 0, u16::MAX]));
    let mut bytes = std::io::Cursor::new(Vec::new());
    image::DynamicImage::ImageRgba16(pixels)
        .write_to(&mut bytes, image::ImageFormat::Png)
        .unwrap();
    let record =
        crate::gel_image::import_gel_image_bytes("alpha", "sixteen.png", bytes.get_ref(), None)
            .unwrap();
    let gray = original_pixels(&record).unwrap();
    assert_eq!(gray.get_pixel(0, 0).0[0], 1.0);
    assert_eq!(gray.get_pixel(1, 1).0[0], 0.0);
    let preview = image::load_from_memory(&STANDARD.decode(&record.preview_png_base64).unwrap())
        .unwrap()
        .to_rgb8();
    assert_eq!(preview.get_pixel(0, 0).0, [255; 3]);
    // Preview coordinates can be rescaled; measurement coordinates cannot.
    assert_eq!(
        preview
            .get_pixel(preview.width() - 1, preview.height() - 1)
            .0,
        [0; 3]
    );
    let exported = crate::gel_image::gel_image_export_preview(&record).unwrap();
    assert_eq!(exported, record.preview_png_base64);
    assert_eq!(
        STANDARD.decode(record.original_base64).unwrap(),
        *bytes.get_ref()
    );
}

#[test]
fn stored_rgb_luminance_preserves_white_and_uses_explicit_display_weights() {
    let rgb = image::RgbImage::from_raw(
        3,
        2,
        vec![
            255, 255, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 128, 128, 128,
        ],
    )
    .unwrap();
    let mut bytes = std::io::Cursor::new(Vec::new());
    image::DynamicImage::ImageRgb8(rgb)
        .write_to(&mut bytes, image::ImageFormat::Png)
        .unwrap();
    let image =
        crate::gel_image::import_gel_image_bytes("rgb", "stored.png", bytes.get_ref(), None)
            .unwrap();
    let gray = original_pixels(&image).unwrap();
    for (pixel, expected) in gray
        .pixels()
        .zip([1.0, 0.0, 0.2126, 0.7152, 0.0722, 128.0 / 255.0])
    {
        assert!((f64::from(pixel.0[0]) - expected).abs() < 1e-7);
        assert!((0.0..=1.0).contains(&pixel.0[0]));
    }
    assert_eq!(
        STANDARD.decode(image.original_base64).unwrap(),
        *bytes.get_ref()
    );
}

#[test]
fn opaque_preview_keeps_legacy_bytes_and_does_not_add_alpha_warnings() {
    let pixels = image::RgbaImage::from_pixel(2, 2, image::Rgba([23, 127, 210, 255]));
    let original = image::DynamicImage::ImageRgba8(pixels);
    let mut bytes = std::io::Cursor::new(Vec::new());
    original
        .write_to(&mut bytes, image::ImageFormat::Png)
        .unwrap();
    let record =
        crate::gel_image::import_gel_image_bytes("opaque", "opaque.png", bytes.get_ref(), None)
            .unwrap();
    let mut legacy = std::io::Cursor::new(Vec::new());
    original
        .thumbnail(1600, 1600)
        .to_rgb8()
        .write_to(&mut legacy, image::ImageFormat::Png)
        .unwrap();
    assert_eq!(
        STANDARD.decode(record.preview_png_base64).unwrap(),
        *legacy.get_ref()
    );
    assert!(
        !record
            .descriptor
            .warnings
            .iter()
            .any(|w| w.contains("alpha_white_v1"))
    );
}

#[test]
fn non_display_float_channels_are_rejected_instead_of_clipped() {
    let rgb = image::ImageBuffer::from_pixel(2, 2, image::Rgb([2.0f32, 0.0, 0.0]));
    let mut bytes = std::io::Cursor::new(Vec::new());
    image::DynamicImage::ImageRgb32F(rgb)
        .write_to(&mut bytes, image::ImageFormat::Tiff)
        .unwrap();
    let image = crate::gel_image::import_gel_image_bytes(
        "hdr",
        "non-display.tiff",
        bytes.get_ref(),
        Some(0),
    )
    .unwrap();
    let error = original_pixels(&image).unwrap_err();
    assert!(error.message.contains("display channels"));
}

#[test]
fn csv_backed_export_checks_original_bytes_and_retains_safe_source_summary() {
    use crate::engine::{Engine, GentleEngine, Operation};
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("source.csv");
    std::fs::write(&path, "left_edge,right_edge,weight\n100,300,2\n").unwrap();
    let request = GelProfileFileRequest {
        path: path.display().to_string(),
        profile_id: "a".into(),
        label: "Source".into(),
        run_id: "run".into(),
        lane_id: "lane".into(),
        size_kind: GelSizeKind::LinearDnaBp,
        signal_units: "test_bin_area".into(),
        measure: "synthetic_declared_integrals".into(),
    };
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::GelWorkbench {
            request: Box::new(GelWorkbenchRequest::ImportProfileFile(request)),
        })
        .unwrap();
    engine
        .apply(Operation::GelWorkbench {
            request: Box::new(GelWorkbenchRequest::Compare(GelCompareRequest {
                right_profile: "a".into(),
                ..self::request()
            })),
        })
        .unwrap();
    let out = dir.path().join("verified.json");
    let request = GelWorkbenchExportRequest {
        path: out.display().to_string(),
        analysis_ids: vec![],
        comparison_ids: vec!["compare".into()],
        title: "CSV original".into(),
    };
    engine
        .apply(Operation::GelWorkbench {
            request: Box::new(GelWorkbenchRequest::Export(request.clone())),
        })
        .unwrap();
    let exported: serde_json::Value = serde_json::from_slice(&std::fs::read(out).unwrap()).unwrap();
    assert_eq!(exported["numeric_sources"][0]["profile_id"], "a");
    assert!(
        exported["numeric_sources"][0]
            .get("original_base64")
            .is_none()
    );
    let mut state = serde_json::to_value(engine.state()).unwrap();
    state["gel_images"]["workbench"]["numeric_sources"]["a"]["original_base64"] =
        STANDARD.encode(b"corrupted").into();
    let mut tampered = GentleEngine::from_state(serde_json::from_value(state).unwrap());
    let mut request = request;
    let bad = dir.path().join("refused.json");
    request.path = bad.display().to_string();
    assert!(
        tampered
            .apply(Operation::GelWorkbench {
                request: Box::new(GelWorkbenchRequest::Export(request))
            })
            .is_err()
    );
    assert!(!bad.exists());
}

#[test]
fn pdf_comparison_only_export_rejects_corrupted_raster_and_replaced_document() {
    use crate::engine::{Engine, GentleEngine, Operation};
    // Owner-approved public report fixtures; reproduction/provenance live in
    // test_files/fixtures/gel_images/tapestation_26r266/README.md.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("source.pdf");
    std::fs::write(
        &path,
        include_bytes!("../../test_files/fixtures/gel_images/tapestation_26r266/before.pdf"),
    )
    .unwrap();
    let mut engine = GentleEngine::new();
    let result = engine
        .apply(Operation::GelWorkbench {
            request: Box::new(GelWorkbenchRequest::ImportSource(GelSourceImportRequest {
                image_id: "pdf".into(),
                path: path.display().to_string(),
                pdf_object: Some(5),
                allow_native_thumbnail: false,
            })),
        })
        .unwrap();
    let mut p = profile("a", &[100.0, 300.0], &[Some(1.0)]);
    p.evidence_kind = GelEvidenceKind::ImageDisplay;
    p.source_sha256 = result.gel_image.unwrap().sha256;
    engine
        .apply(Operation::GelWorkbench {
            request: Box::new(GelWorkbenchRequest::ImportProfile { profile: p }),
        })
        .unwrap();
    engine
        .apply(Operation::GelWorkbench {
            request: Box::new(GelWorkbenchRequest::Compare(GelCompareRequest {
                right_profile: "a".into(),
                ..request()
            })),
        })
        .unwrap();
    for replace_document in [false, true] {
        let mut state = serde_json::to_value(engine.state()).unwrap();
        if replace_document {
            let changed =
                include_bytes!("../../test_files/fixtures/gel_images/tapestation_26r266/after.pdf");
            let source = &mut state["gel_images"]["workbench"]["sources"]["pdf"];
            source["original_base64"] = STANDARD.encode(changed).into();
            source["source_sha256"] = sha256_hex_bytes(changed).into();
        } else {
            state["gel_images"]["images"]["pdf"]["original_base64"] =
                STANDARD.encode(b"bad raster").into();
        }
        let mut changed = GentleEngine::from_state(serde_json::from_value(state).unwrap());
        let out = dir.path().join(format!("refused-{replace_document}.json"));
        assert!(
            changed
                .apply(Operation::GelWorkbench {
                    request: Box::new(GelWorkbenchRequest::Export(GelWorkbenchExportRequest {
                        path: out.display().to_string(),
                        title: "Source verification".into(),
                        analysis_ids: vec![],
                        comparison_ids: vec!["compare".into()],
                    }))
                })
                .is_err()
        );
        assert!(!out.exists());
    }
}

fn selectable_text_operations(
    document: &lopdf::Document,
    bytes: &[u8],
    resources: &lopdf::Dictionary,
    depth: usize,
) -> usize {
    assert!(depth < 16, "bounded generated Form nesting");
    let content = lopdf::content::Content::decode(bytes).unwrap();
    let mut text_count = 0;
    let mut unicode_font = false;
    for operation in content.operations {
        match operation.operator.as_str() {
            "Tf" => {
                let fonts = resources
                    .get_deref(b"Font", document)
                    .unwrap()
                    .as_dict()
                    .unwrap();
                let font = fonts
                    .get_deref(operation.operands[0].as_name().unwrap(), document)
                    .unwrap()
                    .as_dict()
                    .unwrap();
                unicode_font = !font
                    .get_deref(b"ToUnicode", document)
                    .unwrap()
                    .as_stream()
                    .unwrap()
                    .get_plain_content()
                    .unwrap()
                    .is_empty();
            }
            "Tj" | "TJ" => {
                assert!(
                    unicode_font,
                    "real text must retain its Unicode font mapping"
                );
                text_count += 1;
            }
            "Do" => {
                let xobjects = resources
                    .get_deref(b"XObject", document)
                    .unwrap()
                    .as_dict()
                    .unwrap();
                let object = xobjects
                    .get_deref(operation.operands[0].as_name().unwrap(), document)
                    .unwrap()
                    .as_stream()
                    .unwrap();
                if object.dict.get(b"Subtype").unwrap().as_name().unwrap() == b"Form" {
                    let child_resources = object
                        .dict
                        .get_deref(b"Resources", document)
                        .ok()
                        .and_then(|r| r.as_dict().ok())
                        .unwrap_or(resources);
                    text_count += selectable_text_operations(
                        document,
                        &object.get_plain_content().unwrap(),
                        child_resources,
                        depth + 1,
                    );
                }
            }
            _ => {}
        }
    }
    text_count
}

fn saved_tapestation_pdf() -> Vec<u8> {
    use crate::engine::{Engine, GentleEngine, Operation};

    let dir = tempfile::tempdir().unwrap();
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("test_files/fixtures/gel_images/tapestation_26r266");
    let manifest: serde_json::Value =
        serde_json::from_slice(&std::fs::read(root.join("manifest.json")).unwrap()).unwrap();
    let mut engine = GentleEngine::new();
    for (id, analysis_id, lane_id) in [("before", "before_A1", "B1"), ("after", "after_L1", "C1")] {
        let result = engine
            .apply(Operation::GelWorkbench {
                request: Box::new(GelWorkbenchRequest::ImportSource(GelSourceImportRequest {
                    image_id: id.into(),
                    path: root.join(format!("{id}.pdf")).display().to_string(),
                    pdf_object: Some(5),
                    allow_native_thumbnail: false,
                })),
            })
            .unwrap();
        let case = manifest["cases"]
            .as_array()
            .unwrap()
            .iter()
            .find(|case| case["id"] == analysis_id)
            .unwrap();
        // Reuse the existing reviewed manual assignments, not new peak guesses.
        // The direct PDF importer may encode its PNG differently while retaining
        // identical decoded pixels, so bind to the actual imported image hash.
        let mut analysis = crate::gel_image::tests::tapestation_request(&manifest, case);
        analysis.image_sha256 = result.gel_image.unwrap().sha256;
        engine
            .apply(Operation::AnalyzeGelImage {
                request: Box::new(analysis),
            })
            .unwrap();
        engine
            .apply(Operation::GelWorkbench {
                request: Box::new(GelWorkbenchRequest::ExtractProfile(
                    GelExtractProfileRequest {
                        profile_id: id.into(),
                        analysis_id: analysis_id.into(),
                        lane_id: lane_id.into(),
                        run_id: id.into(),
                        dark_bands: true,
                        background: 0.0,
                        masked_ranges: vec![],
                    },
                )),
            })
            .unwrap();
    }
    engine
        .apply(Operation::GelWorkbench {
            request: Box::new(GelWorkbenchRequest::Compare(GelCompareRequest {
                left_profile: "before".into(),
                right_profile: "after".into(),
                range: [100.0, 5000.0],
                short_threshold: Some(500.0),
                ..request()
            })),
        })
        .unwrap();
    // Exercise the serialized project, not the editor's transient presentation.
    let state = serde_json::to_value(engine.state()).unwrap();
    let mut reopened = GentleEngine::from_state(serde_json::from_value(state.clone()).unwrap());
    let output = dir.path().join("before after.pdf");
    let request = GelWorkbenchExportRequest {
        path: output.display().to_string(),
        title: "26R266 before/after: image-display shapes, not detector traces".into(),
        analysis_ids: vec!["before_A1".into(), "after_L1".into()],
        comparison_ids: vec!["compare".into()],
    };
    let result = reopened
        .apply(Operation::GelWorkbench {
            request: Box::new(GelWorkbenchRequest::Export(request.clone())),
        })
        .unwrap();
    let bytes = std::fs::read(&output).unwrap();
    assert_eq!(
        result.gel_workbench.unwrap().exported_sha256,
        Some(sha256_hex_bytes(&bytes))
    );
    assert_eq!(serde_json::to_value(reopened.state()).unwrap(), state);
    let document = lopdf::Document::load_mem(&bytes).unwrap();
    assert_eq!(document.get_pages().len(), 5); // overview, provenance, two details, QQ/comparison
    for id in document.get_pages().values() {
        let page = document.get_dictionary(*id).unwrap();
        let resources = page
            .get_deref(b"Resources", &document)
            .unwrap()
            .as_dict()
            .unwrap();
        // lopdf's text extractor does not follow the Form XObjects used by the
        // vector compositor. Inspect invoked forms, text and Unicode mappings.
        assert!(
            selectable_text_operations(&document, &document.get_page_content(*id), resources, 0)
                > 0
        );
    }
    assert!(
        reopened
            .apply(Operation::GelWorkbench {
                request: Box::new(GelWorkbenchRequest::Export(request)),
            })
            .is_err()
    );
    assert_eq!(std::fs::read(&output).unwrap(), bytes);
    bytes
}

#[test]
fn saved_before_after_pdf_contains_ordered_overview_details_and_qq_without_mutation() {
    saved_tapestation_pdf();
}

#[test]
#[ignore = "Manual PDF visual review; set GENTLE_GEL_WORKBENCH_TEST_REPORT to a new path"]
fn write_tapestation_workbench_pdf_for_visual_review() {
    use std::io::Write;
    let path = std::env::var_os("GENTLE_GEL_WORKBENCH_TEST_REPORT")
        .expect("explicit report path required");
    let bytes = saved_tapestation_pdf();
    std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .unwrap()
        .write_all(&bytes)
        .unwrap();
}
