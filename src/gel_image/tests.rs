//! Entirely synthetic fixtures: a 160x220 grayscale image with rectangular bands
//! at explicitly chosen pixel centers. Recreated in memory by synthetic_png;
//! no private gel images, vendor photos, downloads, or simulated-gel renderer.
//! Tests cover measured-gel math, import, engine persistence, shell parity and exports.

use super::*;
use crate::engine::{Engine, GentleEngine, Operation, ProjectState};
use crate::engine_shell::{execute_shell_command, parse_shell_tokens};
use std::sync::Arc;

fn synthetic_png() -> Vec<u8> {
    let mut pixels = image::GrayImage::from_pixel(160, 220, image::Luma([12]));
    for (x, y) in [(20, 20), (20, 80), (20, 140), (20, 200), (90, 110)] {
        for px in x - 8..=x + 8 {
            for py in y - 1..=y + 1 {
                pixels.put_pixel(px, py, image::Luma([210]));
            }
        }
    }
    let mut bytes = Cursor::new(Vec::new());
    pixels.write_to(&mut bytes, ImageFormat::Png).unwrap();
    bytes.into_inner()
}

fn fixture() -> (GelImageRecord, GelImageAnalysisRequest) {
    let image = import_gel_image_bytes("gel", "synthetic.png", &synthetic_png(), None).unwrap();
    let request = GelImageAnalysisRequest {
        report_id: "sizing".into(),
        image_id: "gel".into(),
        image_sha256: image.descriptor.sha256.clone(),
        size_kind: GelSizeKind::LinearDnaBp,
        migration: GelMigrationDirection::Down,
        lanes: vec![
            GelImageLane {
                id: "marker".into(),
                label: "Custom synthetic ladder".into(),
                min: GelImagePoint { x: 0.0, y: 0.0 },
                max: GelImagePoint { x: 40.0, y: 219.0 },
            },
            GelImageLane {
                id: "sample".into(),
                label: "PCR".into(),
                min: GelImagePoint { x: 60.0, y: 0.0 },
                max: GelImagePoint { x: 130.0, y: 219.0 },
            },
        ],
        ladder: GelImageLadder {
            lane_id: "marker".into(),
            label: "Explicit geometric ladder".into(),
            source: "synthetic: sizes halve every 60 pixels".into(),
            gel_system: None,
            prestained: false,
            bands: [
                (20.0, 2000.0),
                (80.0, 1000.0),
                (140.0, 500.0),
                (200.0, 250.0),
            ]
            .iter()
            .enumerate()
            .map(|(i, &(y, size))| GelLadderBand {
                id: format!("ref{i}"),
                center: GelImagePoint { x: 20.0, y },
                size,
            })
            .collect(),
        },
        sample_bands: vec![GelSampleBand {
            id: "band1".into(),
            lane_id: "sample".into(),
            label: "Unknown PCR band".into(),
            center: GelImagePoint { x: 90.0, y: 110.0 },
            position_half_width_px: Some(30.0),
        }],
    };
    (image, request)
}

#[test]
fn geometric_midpoint_is_not_arithmetic_size_and_localization_is_explicit() {
    let (image, request) = fixture();
    let report = analyze_gel_image(&image.descriptor, &request).unwrap();
    let band = &report.estimates[0];
    assert!((band.estimated_size.unwrap() - (500.0_f64 * 1000.0).sqrt()).abs() < 1e-9);
    let bounds = band.localization_size_bounds.unwrap();
    assert!((bounds[0] - 500.0).abs() < 1e-9 && (bounds[1] - 1000.0).abs() < 1e-9);
    assert_eq!(band.reference_band_ids, ["ref1", "ref2"]);
    assert_eq!(report.algorithm, GEL_CALIBRATION_METHOD);
}

#[test]
fn all_four_migration_directions_preserve_sizes() {
    for direction in [
        GelMigrationDirection::Down,
        GelMigrationDirection::Up,
        GelMigrationDirection::Right,
        GelMigrationDirection::Left,
    ] {
        let (mut image, mut request) = fixture();
        let transform = |p: GelImagePoint| match direction {
            GelMigrationDirection::Down => p,
            GelMigrationDirection::Up => GelImagePoint {
                x: p.x,
                y: 219.0 - p.y,
            },
            GelMigrationDirection::Right => GelImagePoint { x: p.y, y: p.x },
            GelMigrationDirection::Left => GelImagePoint {
                x: 219.0 - p.y,
                y: p.x,
            },
        };
        if matches!(
            direction,
            GelMigrationDirection::Left | GelMigrationDirection::Right
        ) {
            image.descriptor.width = 220;
            image.descriptor.height = 160;
        }
        request.migration = direction;
        for lane in &mut request.lanes {
            let a = transform(lane.min);
            let b = transform(lane.max);
            lane.min = GelImagePoint {
                x: a.x.min(b.x),
                y: a.y.min(b.y),
            };
            lane.max = GelImagePoint {
                x: a.x.max(b.x),
                y: a.y.max(b.y),
            };
        }
        for band in &mut request.ladder.bands {
            band.center = transform(band.center);
        }
        for band in &mut request.sample_bands {
            band.center = transform(band.center);
        }
        let report = analyze_gel_image(&image.descriptor, &request).unwrap();
        assert!((report.estimates[0].estimated_size.unwrap() - 707.1067811865476).abs() < 1e-9);
    }
}

#[test]
fn protein_kda_and_prestained_system_do_not_become_dna_bp() {
    let (image, mut request) = fixture();
    request.size_kind = GelSizeKind::SdsProteinKda;
    for band in &mut request.ladder.bands {
        band.size /= 20.0;
    }
    request.ladder.prestained = true;
    assert!(analyze_gel_image(&image.descriptor, &request).is_err());
    request.ladder.gel_system = Some("synthetic Tris-glycine SDS fixture".into());
    let report = analyze_gel_image(&image.descriptor, &request).unwrap();
    assert!((report.estimates[0].estimated_size.unwrap() - 35.35533905932738).abs() < 1e-9);
    assert!(report.warnings.iter().any(|w| w.contains("Apparent")));
    assert!(
        String::from_utf8(gel_image_analysis_tsv(&report).unwrap())
            .unwrap()
            .contains("sds_protein_kda\tkDa")
    );
}

#[test]
fn outside_range_is_not_extrapolated_and_wide_band_bounds_are_not_invented() {
    let (image, mut request) = fixture();
    request.sample_bands[0].center.y = 215.0;
    let report = analyze_gel_image(&image.descriptor, &request).unwrap();
    assert_eq!(
        report.estimates[0].status,
        GelBandSizingStatus::OutsideCalibratedRange
    );
    assert_eq!(report.estimates[0].estimated_size, None);
    assert_eq!(report.estimates[0].localization_size_bounds, None);
    request.sample_bands[0].center.y = 25.0;
    let report = analyze_gel_image(&image.descriptor, &request).unwrap();
    assert!(report.estimates[0].estimated_size.is_some());
    assert!(report.estimates[0].localization_size_bounds.is_none());
}

#[test]
fn missing_reference_is_not_silently_reassigned_and_two_points_warn() {
    let (image, mut request) = fixture();
    request.ladder.bands.remove(1);
    request.ladder.bands.remove(2);
    request.ladder.bands.reverse();
    let report = analyze_gel_image(&image.descriptor, &request).unwrap();
    assert!((report.estimates[0].estimated_size.unwrap() - 707.1067811865476).abs() < 1e-9);
    assert_eq!(report.estimates[0].reference_band_ids, ["ref0", "ref2"]);
    assert!(
        report
            .warnings
            .iter()
            .any(|warning| warning.contains("Only two"))
    );
}

#[test]
fn invalid_assignments_fail_closed() {
    let (image, request) = fixture();
    let invalid: Vec<Box<dyn Fn(&mut GelImageAnalysisRequest)>> = vec![
        Box::new(|r| r.image_sha256 = "wrong".into()),
        Box::new(|r| r.image_id = "other".into()),
        Box::new(|r| r.ladder.bands.truncate(1)),
        Box::new(|r| r.ladder.bands[1].center.y = 20.0),
        Box::new(|r| r.ladder.bands[1].size = 2000.0),
        Box::new(|r| r.ladder.bands[1].size = -1.0),
        Box::new(|r| r.ladder.bands[1].size = f64::NAN),
        Box::new(|r| r.ladder.bands[0].center.x = 100.0),
        Box::new(|r| r.lanes[1].id = "marker".into()),
        Box::new(|r| r.lanes[1].min.x = -1.0),
        Box::new(|r| r.sample_bands[0].id = "ref0".into()),
        Box::new(|r| r.sample_bands[0].center.y = f64::INFINITY),
        Box::new(|r| r.sample_bands[0].center.x = 159.0),
        Box::new(|r| r.sample_bands[0].lane_id = "marker".into()),
        Box::new(|r| r.sample_bands[0].position_half_width_px = Some(-1.0)),
    ];
    for (i, mutate) in invalid.iter().enumerate() {
        let mut changed = request.clone();
        mutate(&mut changed);
        assert!(
            analyze_gel_image(&image.descriptor, &changed).is_err(),
            "case {i}"
        );
    }
}

#[test]
fn original_16bit_data_is_preserved_and_tiff_page_is_explicit() {
    for format in [ImageFormat::Png, ImageFormat::Tiff] {
        let pixels =
            image::ImageBuffer::from_fn(16, 20, |x, y| image::Luma([(x * 3000 + y * 700) as u16]));
        let mut bytes = Cursor::new(Vec::new());
        image::DynamicImage::ImageLuma16(pixels)
            .write_to(&mut bytes, format)
            .unwrap();
        let page = (format == ImageFormat::Tiff).then_some(0);
        if page.is_some() {
            assert!(import_gel_image_bytes("g", "test.tiff", bytes.get_ref(), None).is_err());
            assert!(import_gel_image_bytes("g", "test.tiff", bytes.get_ref(), Some(1)).is_err());
        }
        let image =
            import_gel_image_bytes("g", "/private/patient/test.image", bytes.get_ref(), page)
                .unwrap();
        assert_eq!(
            STANDARD.decode(&image.original_base64).unwrap(),
            *bytes.get_ref()
        );
        assert_eq!(image.descriptor.color_type, "L16");
        assert_eq!(image.descriptor.source_name, "test.image");
        validate_gel_image_record(&image).unwrap();
    }
}

#[test]
fn malformed_or_tampered_images_are_rejected() {
    assert!(import_gel_image_bytes("g", "g.png", b"not a PNG", None).is_err());
    let (mut image, _) = fixture();
    image.descriptor.sha256 = "0".repeat(64);
    assert!(validate_gel_image_record(&image).is_err());
    let (mut image, _) = fixture();
    image.descriptor.width += 1;
    assert!(validate_gel_image_record(&image).is_err());
}

#[test]
fn export_preview_is_rebuilt_from_original_not_mutable_cached_pixels() {
    let (mut image, _) = fixture();
    let original_preview = image.preview_png_base64.clone();
    image.preview_png_base64 = "not the original picture".into();
    assert_eq!(gel_image_export_preview(&image).unwrap(), original_preview);
}

#[test]
fn jpeg_warns_and_oversized_input_is_rejected_before_decoding() {
    let mut bytes = Cursor::new(Vec::new());
    image::RgbImage::new(4, 4)
        .write_to(&mut bytes, ImageFormat::Jpeg)
        .unwrap();
    let image = import_gel_image_bytes("g", "g.jpg", bytes.get_ref(), None).unwrap();
    assert!(
        image
            .descriptor
            .warnings
            .iter()
            .any(|w| w.contains("lossy"))
    );
    assert!(import_gel_image_bytes("g", "g.png", &vec![0; MAX_GEL_IMAGE_BYTES + 1], None).is_err());
}

#[test]
fn manual_workflow_roundtrips_undo_and_export_without_original_path() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("synthetic.png");
    std::fs::write(&input, synthetic_png()).unwrap();
    let (_, request) = fixture();
    let mut engine = GentleEngine::new();
    let import = Operation::ImportGelImage {
        request: GelImageImportRequest {
            image_id: "gel".into(),
            path: input.to_string_lossy().into(),
            tiff_page: None,
        },
    };
    let result = engine.apply(import.clone()).unwrap();
    assert_eq!(result.gel_image.unwrap().sha256, request.image_sha256);
    assert!(engine.apply(import).is_err());
    let shared = engine.state().clone();
    assert!(Arc::ptr_eq(
        &shared.gel_images.images["gel"],
        &engine.state().gel_images.images["gel"]
    ));
    engine
        .apply(Operation::AnalyzeGelImage {
            request: Box::new(request.clone()),
        })
        .unwrap();
    engine.undo_last_operation().unwrap();
    assert!(engine.state().gel_images.analyses.is_empty());
    engine.redo_last_operation().unwrap();
    let state_path = dir.path().join("project.json");
    engine
        .state()
        .save_to_path(state_path.to_str().unwrap())
        .unwrap();
    std::fs::remove_file(&input).unwrap();
    let mut engine = GentleEngine::from_state(
        ProjectState::load_from_path(state_path.to_str().unwrap()).unwrap(),
    );
    for format in [
        GelImageExportFormat::Json,
        GelImageExportFormat::Tsv,
        GelImageExportFormat::Svg,
    ] {
        let path = dir.path().join(format!("export.{format:?}"));
        let op = Operation::ExportGelImageAnalysis {
            request: GelImageExportRequest {
                report_id: "sizing".into(),
                path: path.to_string_lossy().into(),
                format,
            },
        };
        let undo = engine.undo_available();
        engine.apply(op.clone()).unwrap();
        assert_eq!(engine.undo_available(), undo);
        let exported = std::fs::read_to_string(&path).unwrap();
        assert!(!exported.contains(dir.path().to_str().unwrap()));
        assert!(exported.contains(&request.image_sha256));
        assert!(engine.apply(op).is_err());
        assert_eq!(std::fs::read_to_string(path).unwrap(), exported);
        if format == GelImageExportFormat::Svg {
            let tree =
                resvg::usvg::Tree::from_str(&exported, &resvg::usvg::Options::default()).unwrap();
            assert_eq!(tree.size().width(), 1100.0);
            assert!(exported.contains("data:image/png;base64,"));
            assert!(exported.contains("2000 bp"));
            assert!(exported.contains("707 bp"));
            assert!(exported.contains("Orange: confirmed ladder"));
            assert!(!exported.contains("file://"));
        }
    }
}

#[test]
fn shared_shell_routes_match_typed_engine_and_read_only_inspection() {
    let (image, request) = fixture();
    let mut engine = GentleEngine::new();
    engine
        .state_mut()
        .gel_images
        .images
        .insert("gel".into(), Arc::new(image.clone()));
    let tokens = vec![
        "gel-image".into(),
        "analyze".into(),
        serde_json::to_string(&request).unwrap(),
    ];
    let command = parse_shell_tokens(&tokens).unwrap();
    let result = execute_shell_command(&mut engine, &command).unwrap();
    assert!(result.state_changed);
    let actual: GelImageAnalysisReport =
        serde_json::from_value(result.output["result"]["gel_image_analysis"].clone()).unwrap();
    assert_eq!(
        actual,
        analyze_gel_image(&image.descriptor, &request).unwrap()
    );
    let inspect =
        parse_shell_tokens(&["gel-image".into(), "inspect".into(), "sizing".into()]).unwrap();
    let result = execute_shell_command(&mut engine, &inspect).unwrap();
    assert!(!result.state_changed);
    assert_eq!(
        result.output["result"]["gel_image_analysis"]["request"]["report_id"],
        "sizing"
    );
    let mut invalid = serde_json::to_value(&request).unwrap();
    invalid["automatic_extrapolation"] = true.into();
    assert!(
        parse_shell_tokens(&["gel-image".into(), "analyze".into(), invalid.to_string()]).is_err()
    );
}

#[test]
fn export_rejects_edited_report_and_tracks_external_write_for_safety() {
    let (image, request) = fixture();
    let mut engine = GentleEngine::new();
    engine
        .state_mut()
        .gel_images
        .images
        .insert("gel".into(), Arc::new(image));
    engine
        .apply(Operation::AnalyzeGelImage {
            request: Box::new(request),
        })
        .unwrap();
    Arc::make_mut(
        engine
            .state_mut()
            .gel_images
            .analyses
            .get_mut("sizing")
            .unwrap(),
    )
    .estimates[0]
        .estimated_size = Some(9999.0);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("bad.svg").to_string_lossy().to_string();
    let op = Operation::ExportGelImageAnalysis {
        request: GelImageExportRequest {
            report_id: "sizing".into(),
            path: path.clone(),
            format: GelImageExportFormat::Svg,
        },
    };
    assert_eq!(
        GentleEngine::collect_run_bundle_export_paths(&op),
        vec![path.clone()]
    );
    assert!(engine.apply(op).is_err());
    assert!(!Path::new(&path).exists());
}
