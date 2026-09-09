//! Synthetic manual-editor tests. The in-memory 160x220 image is recreated
//! below; no private image, external standard, or vendor artwork is used.

use super::*;
use crate::engine::{Engine, ProjectState};

fn fixture() -> (Arc<RwLock<GentleEngine>>, GelImageEditor) {
    let mut pixels = image::GrayImage::from_pixel(160, 220, image::Luma([15]));
    for (x, y) in [(20, 20), (20, 80), (20, 140), (20, 200), (90, 110)] {
        for xx in x - 8..=x + 8 {
            for yy in y - 1..=y + 1 {
                pixels.put_pixel(xx, yy, image::Luma([210]));
            }
        }
    }
    let mut bytes = Cursor::new(vec![]);
    pixels
        .write_to(&mut bytes, image::ImageFormat::Png)
        .unwrap();
    let image = Arc::new(
        crate::gel_image::import_gel_image_bytes("gel", "synthetic.png", bytes.get_ref(), None)
            .unwrap(),
    );
    let mut engine = GentleEngine::new();
    engine
        .state_mut()
        .gel_images
        .images
        .insert("gel".into(), image.clone());
    let mut editor = GelImageEditor {
        image_id: "gel".into(),
        image: Some(image.clone()),
        draft: Some(empty_draft(&image.descriptor)),
        ..Default::default()
    };
    editor.add_lane(
        GelImagePoint { x: 0.0, y: 0.0 },
        GelImagePoint { x: 40.0, y: 219.0 },
    );
    editor.add_lane(
        GelImagePoint { x: 60.0, y: 0.0 },
        GelImagePoint { x: 130.0, y: 219.0 },
    );
    editor.draft.as_mut().unwrap().ladder.source =
        "Synthetic reference: size halves every 60 pixels".into();
    editor.tool = MarkTool::Reference;
    for (y, size) in [
        (20.0, 2000.0),
        (80.0, 1000.0),
        (140.0, 500.0),
        (200.0, 250.0),
    ] {
        editor.reference_size = size;
        editor.add_band(GelImagePoint { x: 20.0, y });
    }
    editor.tool = MarkTool::Sample;
    editor.add_band(GelImagePoint { x: 90.0, y: 110.0 });
    (Arc::new(RwLock::new(engine)), editor)
}

#[test]
fn gel_image_editor_maps_pixel_centers_through_zoom_pan_and_resize() {
    for (origin, scale) in [
        (egui::pos2(0.0, 0.0), 0.2),
        (egui::pos2(-420.0, -800.0), 3.0),
        (egui::pos2(90.0, 45.0), 1.0),
    ] {
        let transform = ImageTransform {
            rect: egui::Rect::from_min_size(origin, egui::vec2(160.0, 220.0) * scale),
            width: 160,
            height: 220,
        };
        for p in [
            GelImagePoint { x: 0.0, y: 0.0 },
            GelImagePoint { x: 159.0, y: 219.0 },
            GelImagePoint {
                x: 90.25,
                y: 110.75,
            },
        ] {
            let roundtrip = transform.to_original(transform.to_screen(p)).unwrap();
            assert!((roundtrip.x - p.x).abs() < 0.001 && (roundtrip.y - p.y).abs() < 0.001);
        }
        assert!(
            transform
                .to_original(origin - egui::vec2(1.0, 1.0))
                .is_none()
        );
        assert_eq!(
            transform.to_original(transform.rect.left_top()).unwrap(),
            GelImagePoint { x: 0.0, y: 0.0 }
        );
    }
}

#[test]
fn gel_image_editor_marks_require_explicit_size_and_unambiguous_sample_lane() {
    let (_, mut editor) = fixture();
    editor.tool = MarkTool::Reference;
    let count = editor.draft.as_ref().unwrap().ladder.bands.len();
    editor.add_band(GelImagePoint { x: 20.0, y: 150.0 });
    assert_eq!(editor.draft.as_ref().unwrap().ladder.bands.len(), count);
    editor.reference_size = 125.0;
    editor.add_band(GelImagePoint { x: 90.0, y: 215.0 });
    assert_eq!(editor.draft.as_ref().unwrap().ladder.bands.len(), count);
    editor.add_band(GelImagePoint { x: 20.0, y: 215.0 });
    assert_eq!(editor.reference_size, 0.0);
    editor.tool = MarkTool::Sample;
    editor.add_lane(
        GelImagePoint { x: 65.0, y: 0.0 },
        GelImagePoint { x: 125.0, y: 219.0 },
    );
    editor.selected_lane.clear();
    editor.add_band(GelImagePoint { x: 90.0, y: 110.0 });
    assert_eq!(editor.draft.as_ref().unwrap().sample_bands.len(), 1);
    assert!(editor.status.contains("overlap"));
}

#[test]
fn gel_image_editor_move_stays_in_lane_and_calibrates_via_engine() {
    let (shared, mut editor) = fixture();
    editor.move_selected(GelImagePoint { x: 20.0, y: 100.0 });
    assert_eq!(
        editor.draft.as_ref().unwrap().sample_bands[0].center.x,
        90.0
    );
    editor.move_selected(GelImagePoint { x: 95.0, y: 110.0 });
    let mut request = editor.draft.clone().unwrap();
    request.report_id = "result".into();
    let mut result = execute_job(&shared, GelJob::Analyze(request)).unwrap();
    assert!(
        shared
            .read()
            .unwrap()
            .state()
            .gel_images
            .analyses
            .is_empty()
    );
    shared
        .write()
        .unwrap()
        .commit_detached_execution(result.detached.as_mut().unwrap())
        .unwrap();
    let report = result.report.unwrap();
    assert!((report.estimates[0].estimated_size.unwrap() - 707.1067811865).abs() < 1e-6);
}

#[test]
fn gel_image_editor_drafts_roundtrip_undo_and_do_not_claim_measurements() {
    let (shared, mut editor) = fixture();
    let structural_revision = shared.read().unwrap().structural_revision();
    editor.save_draft(&shared);
    let saved = editor.draft.clone().unwrap();
    {
        let engine = shared.read().unwrap();
        let Operation::SaveGelImageDraft { request } = &engine.operation_log().last().unwrap().op
        else {
            panic!("draft operation")
        };
        assert!(Arc::ptr_eq(
            request,
            &engine.state().gel_images.drafts["gel"]
        ));
        let cloned = engine.operation_log().last().unwrap().op.clone();
        let Operation::SaveGelImageDraft { request: copy } = cloned else {
            panic!("cloned draft")
        };
        assert!(
            Arc::ptr_eq(request, &copy),
            "history must not copy the assignment vectors"
        );
    }
    editor.move_selected(GelImagePoint { x: 95.0, y: 120.0 });
    editor.save_draft(&shared);
    shared.write().unwrap().undo_last_operation().unwrap();
    editor.sync_draft(&shared);
    assert_eq!(editor.draft.as_ref(), Some(&saved));
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("draft.json");
    shared
        .read()
        .unwrap()
        .state()
        .save_to_path(path.to_str().unwrap())
        .unwrap();
    let reopened = ProjectState::load_from_path(path.to_str().unwrap()).unwrap();
    assert_eq!(*reopened.gel_images.drafts["gel"], saved);
    assert!(reopened.gel_images.analyses.is_empty());
    let draft_command = crate::engine_shell::parse_shell_tokens(&[
        "gel-image".into(),
        "save-draft".into(),
        serde_json::to_string(&saved).unwrap(),
    ])
    .unwrap();
    let result =
        crate::engine_shell::execute_shell_command(&mut shared.write().unwrap(), &draft_command)
            .unwrap();
    assert!(result.state_changed);
    assert_eq!(
        shared.read().unwrap().structural_revision(),
        structural_revision,
        "manual gel edits must not invalidate DNA structure caches"
    );
    assert!(result.output["result"]["gel_image_analysis"].is_null());
    let mut bad = saved;
    bad.image_sha256 = "wrong".into();
    assert!(
        shared
            .write()
            .unwrap()
            .apply(Operation::SaveGelImageDraft {
                request: Arc::new(bad)
            })
            .is_err()
    );
}

#[test]
fn gel_image_editor_partial_drafts_are_preserved_but_invalid_sizing_is_rejected() {
    let (shared, mut editor) = fixture();
    editor.draft.as_mut().unwrap().ladder.bands.clear();
    editor.save_draft(&shared);
    assert!(
        shared.read().unwrap().state().gel_images.drafts["gel"]
            .ladder
            .bands
            .is_empty()
    );
    let mut request = editor.draft.clone().unwrap();
    request.report_id = "bad".into();
    assert!(execute_job(&shared, GelJob::Analyze(request.clone())).is_err());
    request.sample_bands[0].center.x = f64::NAN;
    assert!(
        shared
            .write()
            .unwrap()
            .apply(Operation::SaveGelImageDraft {
                request: Arc::new(request)
            })
            .is_err()
    );
}

#[test]
fn gel_image_editor_draft_history_interleaves_with_dna_edits() {
    let (shared, mut editor) = fixture();
    editor.save_draft(&shared);
    let original = editor.draft.clone();
    shared
        .write()
        .unwrap()
        .apply(Operation::CreateSequenceFromText {
            sequence_text: "ACGT".into(),
            output_id: Some("unrelated".into()),
            name: None,
            circular: false,
        })
        .unwrap();
    editor.move_selected(GelImagePoint { x: 90.0, y: 120.0 });
    editor.save_draft(&shared);
    shared.write().unwrap().undo_last_operation().unwrap();
    editor.sync_draft(&shared);
    assert_eq!(editor.draft, original);
    assert!(
        shared
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key("unrelated")
    );
    shared.write().unwrap().undo_last_operation().unwrap();
    assert!(
        !shared
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key("unrelated")
    );
    assert_eq!(
        shared.read().unwrap().state().gel_images.drafts["gel"].as_ref(),
        original.as_ref().unwrap()
    );
    shared.write().unwrap().redo_last_operation().unwrap();
    shared.write().unwrap().redo_last_operation().unwrap();
    editor.sync_draft(&shared);
    assert_eq!(
        editor.draft.as_ref().unwrap().sample_bands[0].center.y,
        120.0
    );
}

fn pending_result(editor: &mut GelImageEditor, result: GelJobResult) {
    let (tx, rx) = mpsc::channel();
    assert!(tx.send(Ok(result)).is_ok());
    editor.job = Some(rx);
}

#[test]
fn gel_image_editor_stale_canceled_and_switched_project_results_never_commit() {
    for case in ["stale", "canceled", "project"] {
        let (shared, mut editor) = fixture();
        let ctx = egui::Context::default();
        editor.poll(&shared, &ctx);
        let mut request = editor.draft.clone().unwrap();
        request.report_id = "pending".into();
        let result = execute_job(&shared, GelJob::Analyze(request)).unwrap();
        pending_result(&mut editor, result);
        match case {
            "stale" => {
                shared
                    .write()
                    .unwrap()
                    .apply(Operation::CreateSequenceFromText {
                        sequence_text: "ACGT".into(),
                        output_id: Some("concurrent".into()),
                        name: None,
                        circular: false,
                    })
                    .unwrap();
                editor.poll(&shared, &ctx);
            }
            "canceled" => {
                editor.canceled = true;
                editor.poll(&shared, &ctx);
            }
            _ => {
                editor.poll(&Arc::new(RwLock::new(GentleEngine::new())), &ctx);
            }
        }
        assert!(
            shared
                .read()
                .unwrap()
                .state()
                .gel_images
                .analyses
                .is_empty(),
            "{case}"
        );
        assert!(editor.report.is_none(), "{case}");
    }
}

fn render_canvas(
    editor: &mut GelImageEditor,
    ctx: &egui::Context,
    events: Vec<egui::Event>,
) -> ImageTransform {
    let mut transform = None;
    let mut output = ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                Pos2::ZERO,
                egui::vec2(900.0, 700.0),
            )),
            events,
            ..Default::default()
        },
        |ui| {
            transform = editor.canvas(ui);
        },
    );
    output.textures_delta.clear();
    transform.unwrap()
}

#[test]
fn gel_image_editor_real_pointer_click_uses_original_image_coordinates() {
    let (_, mut editor) = fixture();
    let ctx = egui::Context::default();
    let (_, pixels) = prepare_preview(editor.image.clone().unwrap()).unwrap();
    editor.texture = Some(ctx.load_texture("test", pixels, egui::TextureOptions::LINEAR));
    editor.tool = MarkTool::Sample;
    let transform = render_canvas(&mut editor, &ctx, vec![]);
    let target = transform.to_screen(GelImagePoint { x: 90.0, y: 170.0 });
    render_canvas(
        &mut editor,
        &ctx,
        vec![
            egui::Event::PointerMoved(target),
            egui::Event::PointerButton {
                pos: target,
                button: egui::PointerButton::Primary,
                pressed: true,
                modifiers: Default::default(),
            },
        ],
    );
    render_canvas(
        &mut editor,
        &ctx,
        vec![egui::Event::PointerButton {
            pos: target,
            button: egui::PointerButton::Primary,
            pressed: false,
            modifiers: Default::default(),
        }],
    );
    let bands = &editor.draft.as_ref().unwrap().sample_bands;
    assert_eq!(bands.len(), 2);
    assert!((bands[1].center.x - 90.0).abs() < 0.001 && (bands[1].center.y - 170.0).abs() < 0.001);
    editor.tool = MarkTool::Move;
    let moved = transform.to_screen(GelImagePoint { x: 95.0, y: 180.0 });
    render_canvas(
        &mut editor,
        &ctx,
        vec![
            egui::Event::PointerMoved(target),
            egui::Event::PointerButton {
                pos: target,
                button: egui::PointerButton::Primary,
                pressed: true,
                modifiers: Default::default(),
            },
        ],
    );
    render_canvas(&mut editor, &ctx, vec![egui::Event::PointerMoved(moved)]);
    render_canvas(
        &mut editor,
        &ctx,
        vec![egui::Event::PointerButton {
            pos: moved,
            button: egui::PointerButton::Primary,
            pressed: false,
            modifiers: Default::default(),
        }],
    );
    let moved_band = &editor.draft.as_ref().unwrap().sample_bands[1];
    assert!(
        (moved_band.center.x - 95.0).abs() < 0.001 && (moved_band.center.y - 180.0).abs() < 0.001
    );
}

#[test]
fn gel_image_editor_idle_frames_reuse_texture_and_do_not_mutate_project() {
    let (shared, mut editor) = fixture();
    let ctx = egui::Context::default();
    let (_, pixels) = prepare_preview(editor.image.clone().unwrap()).unwrap();
    editor.texture = Some(ctx.load_texture("test", pixels, egui::TextureOptions::LINEAR));
    editor.save_draft(&shared);
    let revision = shared.read().unwrap().mutation_revision();
    let texture = editor.texture.as_ref().unwrap().id();
    for frame in 0..3 {
        let mut output = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    Pos2::ZERO,
                    egui::vec2(1200.0, 1000.0),
                )),
                ..Default::default()
            },
            |ui| editor.contents(ui, &shared),
        );
        if frame > 0 {
            assert!(
                output.textures_delta.set.is_empty(),
                "idle frame uploaded a texture"
            );
        }
        output.textures_delta.clear();
        assert_eq!(editor.texture.as_ref().unwrap().id(), texture);
        assert!(editor.job.is_none());
    }
    assert_eq!(shared.read().unwrap().mutation_revision(), revision);
}

#[test]
fn gel_image_editor_menu_shell_and_window_inventory_share_one_target() {
    let mut app = GENtleApp::default();
    let parse = |action: &str| {
        crate::engine_shell::parse_shell_tokens(&[
            "ui".into(),
            action.into(),
            "gel-image-editor".into(),
        ])
        .unwrap()
    };
    app.try_apply_shell_ui_intent(&parse("open")).unwrap();
    assert!(app.gel_image_editor.open);
    assert!(
        app.collect_open_window_entries()
            .iter()
            .any(|e| e.viewport_id == GENtleApp::gel_image_viewport_id())
    );
    app.try_apply_shell_ui_intent(&parse("close")).unwrap();
    assert!(!app.gel_image_editor.open);
    app.try_apply_shell_ui_intent(&parse("focus")).unwrap();
    assert!(app.gel_image_editor.open);
    app.engine = Arc::new(RwLock::new(GentleEngine::new()));
    app.try_apply_shell_ui_intent(&parse("open")).unwrap();
    app.gel_image_editor
        .poll(&app.engine, &egui::Context::default());
    assert!(
        app.gel_image_editor.open,
        "opening immediately after a project switch must not be lost"
    );
    let target = UiIntentTarget::GelImageEditor;
    assert_eq!(target.menu_path(), "Patterns");
    assert!(target.arguments().is_empty());
    assert!(
        crate::engine_shell::parse_shell_tokens(&[
            "ui".into(),
            "open".into(),
            "gel-image-editor".into(),
            "--genome-id".into(),
            "unused".into()
        ])
        .is_err()
    );
}
