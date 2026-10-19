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

fn detection_result(
    shared: &Arc<RwLock<GentleEngine>>,
    editor: &mut GelImageEditor,
) -> GelJobResult {
    // Reuses the synthetic two-lane image above, not a vendor reference.
    editor.workbench.detection_lanes = editor
        .draft
        .as_ref()
        .unwrap()
        .lanes
        .iter()
        .map(|lane| lane.id.clone())
        .collect();
    editor.save_draft(shared);
    let request = editor.detection_request().unwrap();
    assert_eq!(
        request.lanes.len(),
        1,
        "ladder is excluded even if selected"
    );
    assert_ne!(
        request.lanes[0].id,
        editor.draft.as_ref().unwrap().ladder.lane_id
    );
    let result = execute_job(
        shared,
        GelJob::Detect {
            request,
            draft: Box::new(editor.draft.clone().unwrap()),
        },
    )
    .unwrap();
    assert!(
        result.detached.is_none(),
        "detection must not commit a read-only engine"
    );
    assert!(
        !result
            .suggestions
            .as_ref()
            .unwrap()
            .suggestions
            .bands
            .is_empty()
    );
    result
}

#[test]
fn gel_workbench_pending_suggestions_require_exact_draft_and_request_identity() {
    for case in [
        "draft",
        "detector",
        "selection",
        "stored_draft",
        "cancel",
        "project",
    ] {
        let (shared, mut editor) = fixture();
        let ctx = egui::Context::default();
        editor.poll(&shared, &ctx);
        let result = detection_result(&shared, &mut editor);
        let original_count = editor.draft.as_ref().unwrap().sample_bands.len();
        pending_result(&mut editor, result);
        match case {
            "draft" => {
                editor.draft.as_mut().unwrap().sample_bands[0].label = "Edited while pending".into()
            }
            "detector" => editor.workbench.prominence = 0.2,
            "selection" => editor.workbench.detection_lanes.clear(),
            "stored_draft" => {
                let mut draft = editor.draft.clone().unwrap();
                draft.lanes[1].label = "External draft edit".into();
                shared
                    .write()
                    .unwrap()
                    .apply(Operation::SaveGelImageDraft {
                        request: Arc::new(draft),
                    })
                    .unwrap();
            }
            "cancel" => editor.canceled = true,
            _ => {}
        }
        if case == "project" {
            editor.poll(&Arc::new(RwLock::new(GentleEngine::new())), &ctx);
        } else {
            editor.poll(&shared, &ctx);
        }
        assert!(editor.workbench.pending.is_none(), "{case}");
        assert_eq!(
            shared.read().unwrap().state().gel_images.drafts["gel"]
                .sample_bands
                .len(),
            original_count,
            "{case}"
        );
    }
}

#[test]
fn gel_workbench_reviewed_proposals_append_only_explicit_sample_marks() {
    let (shared, mut editor) = fixture();
    let ctx = egui::Context::default();
    editor.poll(&shared, &ctx);
    let result = detection_result(&shared, &mut editor);
    let ladder = editor.draft.as_ref().unwrap().ladder.clone();
    let count = editor.draft.as_ref().unwrap().sample_bands.len();
    let revision = shared.read().unwrap().mutation_revision();
    pending_result(&mut editor, result);
    editor.poll(&shared, &ctx);
    assert_eq!(shared.read().unwrap().mutation_revision(), revision);
    assert_eq!(
        editor.draft.as_ref().unwrap().sample_bands.len(),
        count,
        "polling proposals must not accept them"
    );
    editor.accept_suggestions(&shared);
    assert_eq!(
        editor.draft.as_ref().unwrap().sample_bands.len(),
        count,
        "no selection means no acceptance"
    );
    let pending = editor.workbench.pending.as_ref().unwrap();
    editor
        .workbench
        .selected_suggestions
        .insert(pending.suggestions.bands[0].id.clone());
    editor.accept_suggestions(&shared);
    let draft = editor.draft.as_ref().unwrap();
    assert_eq!(draft.ladder, ladder);
    assert_eq!(draft.sample_bands.len(), count + 1);
    assert!(
        draft
            .sample_bands
            .iter()
            .all(|band| band.lane_id != ladder.lane_id)
    );
    assert_eq!(
        &*shared.read().unwrap().state().gel_images.drafts["gel"],
        draft
    );
    assert!(editor.workbench.pending.is_none());
    assert!(
        shared
            .read()
            .unwrap()
            .state()
            .gel_images
            .analyses
            .is_empty(),
        "accepted drafts are not sizing reports"
    );
}

#[test]
fn gel_workbench_suggestions_cannot_be_accepted_after_review_target_changes() {
    for change_detector in [false, true] {
        let (shared, mut editor) = fixture();
        let ctx = egui::Context::default();
        editor.poll(&shared, &ctx);
        let result = detection_result(&shared, &mut editor);
        pending_result(&mut editor, result);
        editor.poll(&shared, &ctx);
        let pending = editor.workbench.pending.as_ref().unwrap();
        editor
            .workbench
            .selected_suggestions
            .insert(pending.suggestions.bands[0].id.clone());
        let saved = shared.read().unwrap().state().gel_images.drafts["gel"].clone();
        let revision = shared.read().unwrap().mutation_revision();
        if change_detector {
            editor.workbench.dark_bands = !editor.workbench.dark_bands;
        } else {
            editor.draft.as_mut().unwrap().lanes[1].max.x -= 1.0;
        }
        editor.accept_suggestions(&shared);
        assert_eq!(shared.read().unwrap().mutation_revision(), revision);
        assert!(Arc::ptr_eq(
            &saved,
            &shared.read().unwrap().state().gel_images.drafts["gel"]
        ));
        assert!(editor.workbench.pending.is_none());
    }
}

#[test]
fn gel_workbench_empty_inventory_disables_dependent_actions_without_an_image() {
    let inventory = workbench::Inventory::default();
    let mut editor = GelImageEditor::default();
    assert!(editor.detection_request().is_none());
    assert!(editor.workbench.compare_request(&inventory).is_none());
    assert!(editor.workbench.extraction_request(&inventory).is_none());
    assert!(editor.workbench.algebra_request(&inventory).is_none());
    assert!(editor.workbench.export_request(&inventory).is_none());
    editor.workbench.left = "missing".into();
    editor.workbench.right = "also_missing".into();
    editor.workbench.export_analyses.push("missing".into());
    assert!(editor.workbench.compare_request(&inventory).is_none());
    assert!(editor.workbench.export_request(&inventory).is_none());
    let shared = Arc::new(RwLock::new(GentleEngine::new()));
    let revision = shared.read().unwrap().mutation_revision();
    let ctx = egui::Context::default();
    let mut output = ctx.run_ui(egui::RawInput::default(), |ui| editor.contents(ui, &shared));
    output.textures_delta.clear();
    assert!(editor.job.is_none());
    assert_eq!(shared.read().unwrap().mutation_revision(), revision);
}

// Hand-crafted integrated-bin profiles, recreated here. They model no real
// instrument or biology and exercise only GUI request/result plumbing.
fn numeric_profile(id: &str, run: &str, hash_digit: char, weights: [f64; 2]) -> GelProfile {
    GelProfile {
        id: id.into(),
        label: format!("Synthetic {id}"),
        run_id: run.into(),
        lane_id: "lane_2".into(),
        source_sha256: hash_digit.to_string().repeat(64),
        evidence_kind: GelEvidenceKind::VendorTrace,
        size_kind: GelSizeKind::LinearDnaBp,
        signal_units: "synthetic integrated signal".into(),
        edges: vec![100.0, 200.0, 400.0],
        weights: weights.into_iter().map(Some).collect(),
        measure: "synthetic_integrated_bins_v1".into(),
        processing: vec!["Hand-crafted GUI test bins".into()],
        warnings: vec![],
    }
}

fn import_numeric_pair(
    shared: &Arc<RwLock<GentleEngine>>,
    editor: &mut GelImageEditor,
    ctx: &egui::Context,
) {
    for profile in [
        numeric_profile("left", "before_run", 'a', [2.0, 6.0]),
        numeric_profile("right", "after_run", 'b', [4.0, 1.0]),
    ] {
        let result = execute_job(
            shared,
            GelJob::Workbench(GelWorkbenchRequest::ImportProfile { profile }),
        )
        .unwrap();
        assert!(result.detached.is_some());
        pending_result(editor, result);
        editor.poll(shared, ctx);
    }
    editor.workbench.left = "left".into();
    editor.workbench.right = "right".into();
    editor.workbench.range = [100.0, 400.0];
}

#[test]
fn gel_workbench_cross_run_comparison_is_saved_and_idle_frames_do_not_rescore() {
    let shared = Arc::new(RwLock::new(GentleEngine::new()));
    let ctx = egui::Context::default();
    let mut editor = GelImageEditor {
        help: false,
        ..Default::default()
    };
    import_numeric_pair(&shared, &mut editor, &ctx);
    assert!(
        editor.image.is_none(),
        "numeric-only projects must not need an image texture"
    );
    let inventory = workbench::Inventory {
        analyses: BTreeMap::new(),
        store: shared.read().unwrap().state().gel_images.workbench.clone(),
    };
    let request = editor.workbench.compare_request(&inventory).unwrap();
    let result = execute_job(&shared, GelJob::Workbench(request)).unwrap();
    assert!(
        shared
            .read()
            .unwrap()
            .state()
            .gel_images
            .workbench
            .comparisons
            .is_empty()
    );
    pending_result(&mut editor, result);
    editor.poll(&shared, &ctx);
    let report = shared
        .read()
        .unwrap()
        .state()
        .gel_images
        .workbench
        .comparisons["comparison_1"]
        .clone();
    assert_eq!(report.left.run_id, "before_run");
    assert_eq!(report.right.run_id, "after_run");
    assert_ne!(report.left.source_sha256, report.right.source_sha256);
    assert!(report.p_value.is_none());
    let revision = shared.read().unwrap().mutation_revision();
    let history = shared.read().unwrap().operation_log().len();
    editor.workbench.range = [1.0, 2.0];
    editor.workbench.left.clear();
    for _ in 0..3 {
        let mut output = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    Pos2::ZERO,
                    egui::vec2(1200.0, 1400.0),
                )),
                ..Default::default()
            },
            |ui| editor.contents(ui, &shared),
        );
        output.textures_delta.clear();
        assert!(editor.job.is_none());
    }
    assert!(Arc::ptr_eq(
        &report,
        &shared
            .read()
            .unwrap()
            .state()
            .gel_images
            .workbench
            .comparisons["comparison_1"]
    ));
    assert_eq!(shared.read().unwrap().mutation_revision(), revision);
    assert_eq!(shared.read().unwrap().operation_log().len(), history);
    assert_eq!(report.request.range, [100.0, 400.0]);
}

#[test]
fn gel_workbench_signed_algebra_uses_explicit_unclipped_coefficients() {
    let shared = Arc::new(RwLock::new(GentleEngine::new()));
    let ctx = egui::Context::default();
    let mut editor = GelImageEditor::default();
    import_numeric_pair(&shared, &mut editor, &ctx);
    editor.workbench.coefficients = [1.0, -2.0];
    let inventory = workbench::Inventory {
        analyses: BTreeMap::new(),
        store: shared.read().unwrap().state().gel_images.workbench.clone(),
    };
    let request = editor.workbench.algebra_request(&inventory).unwrap();
    let GelWorkbenchRequest::Algebra(algebra) = &request else {
        panic!("typed algebra");
    };
    assert_eq!(algebra.terms[0].coefficient, 1.0);
    assert_eq!(algebra.terms[1].coefficient, -2.0);
    let result = execute_job(&shared, GelJob::Workbench(request)).unwrap();
    pending_result(&mut editor, result);
    editor.poll(&shared, &ctx);
    let store = shared.read().unwrap().state().gel_images.workbench.clone();
    let output = &store.profiles["residual_1"];
    assert_eq!(output.weights, vec![Some(-6.0), Some(4.0)]);
    assert_eq!(output.evidence_kind, GelEvidenceKind::SignedResidual);
    editor.workbench.left = output.id.clone();
    assert!(
        editor
            .workbench
            .compare_request(&workbench::Inventory {
                analyses: BTreeMap::new(),
                store
            })
            .is_none()
    );
}

fn algebra_editor() -> (Arc<RwLock<GentleEngine>>, GelImageEditor, egui::Context) {
    let shared = Arc::new(RwLock::new(GentleEngine::new()));
    let ctx = egui::Context::default();
    let mut editor = GelImageEditor {
        open: true,
        help: false,
        ..Default::default()
    };
    import_numeric_pair(&shared, &mut editor, &ctx);
    let inventory = workbench::Inventory {
        store: shared.read().unwrap().state().gel_images.workbench.clone(),
        ..Default::default()
    };
    let result = execute_job(
        &shared,
        GelJob::Workbench(editor.workbench.algebra_request(&inventory).unwrap()),
    )
    .unwrap();
    pending_result(&mut editor, result);
    editor.poll(&shared, &ctx);
    assert_eq!(editor.workbench.selected_algebra, "residual_1");
    (shared, editor, ctx)
}

#[test]
fn gel_workbench_algebra_inspection_rasterizes_off_frame_and_reuses_shared_svg_result() {
    let (shared, mut editor, ctx) = algebra_editor();
    let result = execute_job(
        &shared,
        GelJob::Workbench(GelWorkbenchRequest::InspectAlgebra {
            profile_id: "residual_1".into(),
        }),
    )
    .unwrap();
    assert!(result.detached.is_none());
    let cached = result.algebra_inspection.as_ref().unwrap();
    assert!(cached.pixels.is_some());
    assert!(cached.texture.is_none());
    assert_eq!(
        &*cached.svg,
        gentle_render::gel_workbench::render_algebra_svg(&cached.presentation).as_bytes()
    );
    let svg = cached.svg.clone();
    pending_result(&mut editor, result);
    editor.poll(&shared, &ctx);
    assert!(
        editor
            .workbench
            .algebra_inspection
            .as_ref()
            .unwrap()
            .pixels
            .is_none()
    );
    assert!(
        editor
            .workbench
            .algebra_inspection
            .as_ref()
            .unwrap()
            .texture
            .is_some()
    );
    let revision = shared.read().unwrap().mutation_revision();
    let history = shared.read().unwrap().operation_log().len();
    for _ in 0..3 {
        let mut output = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    Pos2::ZERO,
                    egui::vec2(1200.0, 1400.0),
                )),
                ..Default::default()
            },
            |ui| editor.contents(ui, &shared),
        );
        output.textures_delta.clear();
        assert!(editor.job.is_none());
        assert!(Arc::ptr_eq(
            &svg,
            &editor.workbench.algebra_inspection.as_ref().unwrap().svg
        ));
    }
    assert_eq!(shared.read().unwrap().mutation_revision(), revision);
    assert_eq!(shared.read().unwrap().operation_log().len(), history);
    let inventory = workbench::Inventory {
        store: shared.read().unwrap().state().gel_images.workbench.clone(),
        ..Default::default()
    };
    editor.workbench.export_algebra = vec!["residual_1".into()];
    let GelWorkbenchRequest::Export(request) = editor.workbench.export_request(&inventory).unwrap()
    else {
        panic!("export");
    };
    assert_eq!(request.algebra_ids, vec!["residual_1"]);
}

#[test]
fn gel_workbench_algebra_inspection_rejects_stale_canceled_closed_and_changed_selection() {
    for mode in 0..5 {
        let (shared, mut editor, ctx) = algebra_editor();
        let result = execute_job(
            &shared,
            GelJob::Workbench(GelWorkbenchRequest::InspectAlgebra {
                profile_id: "residual_1".into(),
            }),
        )
        .unwrap();
        pending_result(&mut editor, result);
        match mode {
            0 => editor.canceled = true,
            1 => editor.open = false,
            2 => editor.workbench.selected_algebra = "other".into(),
            3 => {
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
            }
            _ => {
                let other = Arc::new(RwLock::new(GentleEngine::new()));
                editor.poll(&other, &ctx);
                assert!(editor.workbench.algebra_inspection.is_none());
                continue;
            }
        }
        editor.poll(&shared, &ctx);
        assert!(editor.workbench.algebra_inspection.is_none());
    }
}

#[test]
fn gel_workbench_extraction_binds_saved_calibration_lane_run_and_background() {
    let (shared, mut editor) = fixture();
    let ctx = egui::Context::default();
    editor.poll(&shared, &ctx);
    let mut analysis = editor.draft.clone().unwrap();
    analysis.report_id = "calibration".into();
    let result = execute_job(&shared, GelJob::Analyze(analysis)).unwrap();
    pending_result(&mut editor, result);
    editor.poll(&shared, &ctx);
    let report = shared.read().unwrap().state().gel_images.analyses["calibration"].clone();
    let inventory = workbench::Inventory {
        analyses: [("calibration".into(), report.clone())]
            .into_iter()
            .collect(),
        ..Default::default()
    };
    editor.workbench.analysis_id = "calibration".into();
    editor.workbench.lane_id = report.request.ladder.lane_id.clone();
    editor.workbench.run_id = "Explicit run label".into();
    editor.workbench.background = 0.1;
    assert!(
        editor.workbench.extraction_request(&inventory).is_none(),
        "ladder is not a sample profile"
    );
    editor.workbench.lane_id = report.request.lanes[1].id.clone();
    editor.draft.as_mut().unwrap().migration = GelMigrationDirection::Left;
    let request = editor.workbench.extraction_request(&inventory).unwrap();
    let GelWorkbenchRequest::ExtractProfile(extract) = &request else {
        panic!("typed extraction");
    };
    assert_eq!(extract.analysis_id, "calibration");
    assert_eq!(extract.run_id, "Explicit run label");
    assert!(!extract.dark_bands);
    assert_eq!(extract.background, 0.1);
    assert!(extract.masked_ranges.is_empty());
    let result = execute_job(&shared, GelJob::Workbench(request)).unwrap();
    assert!(result.detached.is_some());
    let profile = result.workbench.as_ref().unwrap().profile.as_ref().unwrap();
    assert_eq!(profile.evidence_kind, GelEvidenceKind::ImageDisplay);
    assert_eq!(profile.lane_id, report.request.lanes[1].id);
    assert!(profile.weights.iter().flatten().any(|w| *w > 0.0));
    pending_result(&mut editor, result);
    editor.poll(&shared, &ctx);
    assert!(
        shared
            .read()
            .unwrap()
            .state()
            .gel_images
            .workbench
            .profiles
            .contains_key("profile_1")
    );
}

#[test]
fn gel_workbench_mutating_jobs_reject_stale_and_canceled_results() {
    for canceled in [false, true] {
        let shared = Arc::new(RwLock::new(GentleEngine::new()));
        let ctx = egui::Context::default();
        let mut editor = GelImageEditor::default();
        editor.poll(&shared, &ctx);
        let result = execute_job(
            &shared,
            GelJob::Workbench(GelWorkbenchRequest::ImportProfile {
                profile: numeric_profile("pending", "synthetic run", 'a', [2.0, 6.0]),
            }),
        )
        .unwrap();
        pending_result(&mut editor, result);
        if canceled {
            editor.canceled = true;
        } else {
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
        }
        editor.poll(&shared, &ctx);
        assert!(
            shared
                .read()
                .unwrap()
                .state()
                .gel_images
                .workbench
                .profiles
                .is_empty()
        );
    }
}

#[test]
fn gel_workbench_export_is_external_and_preserves_explicit_analysis_order() {
    let (shared, editor) = fixture();
    let mut analysis = editor.draft.clone().unwrap();
    analysis.report_id = "a_before".into();
    let first = execute_job(&shared, GelJob::Analyze(analysis.clone()))
        .unwrap()
        .report
        .unwrap();
    analysis.report_id = "z_after".into();
    let second = execute_job(&shared, GelJob::Analyze(analysis))
        .unwrap()
        .report
        .unwrap();
    let inventory = workbench::Inventory {
        analyses: [("a_before".into(), first), ("z_after".into(), second)]
            .into_iter()
            .collect(),
        ..Default::default()
    };
    let mut state = WorkbenchUi::default();
    state.export_analyses = vec!["z_after".into(), "a_before".into()];
    let request = state.export_request(&inventory).unwrap();
    let GelWorkbenchRequest::Export(export) = &request else {
        panic!("typed PDF export");
    };
    assert_eq!(export.analysis_ids, ["z_after", "a_before"]);
    assert!(!request.mutates());
    assert!(GelJob::Workbench(request).is_external_export());
    assert!(
        GelJob::Export(GelImageExportRequest {
            report_id: "a_before".into(),
            path: "unused.json".into(),
            format: GelImageExportFormat::Json,
        })
        .is_external_export()
    );
    assert!(!GelJob::Workbench(GelWorkbenchRequest::Inspect {}).is_external_export());
    state.export_analyses.push("removed_report".into());
    assert!(state.export_request(&inventory).is_none());
}

#[test]
fn gel_workbench_export_job_tracking_cannot_offer_cancellation_or_commit() {
    let shared = Arc::new(RwLock::new(GentleEngine::new()));
    let revision = shared.read().unwrap().mutation_revision();
    let ctx = egui::Context::default();
    let mut editor = GelImageEditor::default();
    // Invalid suffix stops before rendering or writing; this tests dispatch
    // classification without requiring a native PDF renderer or a file dialog.
    editor.start(
        &shared,
        &ctx,
        GelJob::Workbench(GelWorkbenchRequest::Export(GelWorkbenchExportRequest {
            path: "not-a-pdf.txt".into(),
            analysis_ids: vec![],
            comparison_ids: vec![],
            title: "Synthetic export".into(),
            algebra_ids: vec![],
        })),
    );
    assert!(editor.export_running);
    assert!(editor.job.is_some());
    assert!(!editor.canceled);
    assert_eq!(shared.read().unwrap().mutation_revision(), revision);
    pending_result(
        &mut editor,
        GelJobResult {
            detached: None,
            preview: None,
            report: None,
            source: None,
            suggestions: None,
            workbench: Some(Box::new(GelWorkbenchResult {
                exported_sha256: Some("a".repeat(64)),
                ..Default::default()
            })),
            algebra_inspection: None,
            message: "Synthetic external-export receipt; no file produced".into(),
        },
    );
    editor.poll(&shared, &ctx);
    assert!(!editor.export_running);
    assert!(editor.status.contains("no project commit"));
    assert_eq!(shared.read().unwrap().mutation_revision(), revision);
}

#[test]
fn gel_workbench_qq_coordinates_use_saved_size_values_and_upward_right_axis() {
    let rect = egui::Rect::from_min_size(egui::pos2(10.0, 20.0), egui::vec2(200.0, 200.0));
    let range = [100.0, 500.0];
    assert_eq!(
        workbench::qq_position(rect, range, 100.0, 100.0),
        Some(rect.left_bottom())
    );
    assert_eq!(
        workbench::qq_position(rect, range, 500.0, 500.0),
        Some(rect.right_top())
    );
    assert_eq!(
        workbench::qq_position(rect, range, 200.0, 400.0),
        Some(egui::pos2(60.0, 70.0))
    );
    assert_eq!(
        workbench::qq_position(rect, range, 300.0, 300.0),
        Some(rect.center())
    );
    for (bounds, left, right) in [
        (range, f64::NAN, 200.0),
        (range, 200.0, 501.0),
        (range, 99.0, 200.0),
        ([500.0, 100.0], 200.0, 300.0),
    ] {
        assert!(
            workbench::qq_position(rect, bounds, left, right).is_none(),
            "no clipping into the plot"
        );
    }
}

#[test]
fn gel_workbench_native_import_requires_pdf_object_and_thumbnail_opt_in() {
    let mut request = GelWorkbenchRequest::ImportSource(GelSourceImportRequest {
        image_id: "source".into(),
        path: "synthetic.pdf".into(),
        pdf_object: None,
        allow_native_thumbnail: false,
    });
    assert!(
        workbench::prepare_file_request(&mut request)
            .unwrap_err()
            .contains("explicit PDF")
    );
    let GelWorkbenchRequest::ImportSource(source) = &mut request else {
        unreachable!()
    };
    source.pdf_object = Some(17);
    assert!(workbench::prepare_file_request(&mut request).is_ok());
    let GelWorkbenchRequest::ImportSource(source) = &mut request else {
        unreachable!()
    };
    source.path = "synthetic.D5000".into();
    assert!(
        workbench::prepare_file_request(&mut request)
            .unwrap_err()
            .contains("Encrypted numerical")
    );
    let GelWorkbenchRequest::ImportSource(source) = &mut request else {
        unreachable!()
    };
    source.allow_native_thumbnail = true;
    assert!(workbench::prepare_file_request(&mut request).is_ok());
}

#[test]
fn gel_workbench_numeric_file_request_roundtrips_without_guessing_integration() {
    let request = GelWorkbenchRequest::ImportProfileFile(GelProfileFileRequest {
        path: "synthetic bins.csv".into(),
        profile_id: "bins".into(),
        label: "Synthetic bins".into(),
        run_id: "run_a".into(),
        lane_id: "lane_2".into(),
        size_kind: GelSizeKind::LinearDnaBp,
        signal_units: "synthetic signal".into(),
        measure: "caller_declared_integrated_bins".into(),
    });
    let json = serde_json::to_string(&request).unwrap();
    assert!(json.contains("import_profile_file"));
    let mut parsed = serde_json::from_str::<GelWorkbenchRequest>(&json).unwrap();
    assert_eq!(request, parsed);
    assert!(parsed.mutates());
    assert!(workbench::prepare_file_request(&mut parsed).is_ok());
}

#[test]
fn gel_workbench_removed_image_preview_is_not_adopted() {
    let (shared, mut editor) = fixture();
    let ctx = egui::Context::default();
    editor.poll(&shared, &ctx);
    let result = execute_job(&shared, GelJob::Preview(editor.image.clone().unwrap())).unwrap();
    shared
        .write()
        .unwrap()
        .state_mut()
        .gel_images
        .images
        .remove("gel");
    pending_result(&mut editor, result);
    editor.poll(&shared, &ctx);
    assert!(editor.image.is_none());
    assert!(editor.texture.is_none());
    assert!(editor.status.contains("stale preview"));
}
