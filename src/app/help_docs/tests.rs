//! Lazy Help regressions using the public catalog and temporary, hand-written
//! Markdown/SVG fixtures recreated below. No user files or font settings change.

use super::*;
use crate::gui_profiler::startup_trace;

#[test]
fn gui_manual_pages_reuse_index_and_render_cache_until_source_changes() {
    let mut app = GENtleApp::default();
    app.help_gui_markdown =
        "# Manual\n\nOverview\n\n## First\n\nFIRST_BODY\n\n## Second\n\nSECOND_BODY\n".into();
    let source = app.help_gui_markdown.clone();
    let index = app.ensure_gui_manual_pages();
    app.select_gui_manual_page(1);
    let first = app.rendered_help_markdown_for(HelpDoc::Gui);
    let second = app.rendered_help_markdown_for(HelpDoc::Gui);
    assert!(Arc::ptr_eq(&first, &second));
    assert!(Arc::ptr_eq(&index, &app.ensure_gui_manual_pages()));
    assert!(first.contains("FIRST_BODY"));
    assert!(!first.contains("SECOND_BODY"));
    app.select_gui_manual_page(2);
    assert!(
        app.rendered_help_markdown_for(HelpDoc::Gui)
            .contains("SECOND_BODY")
    );
    assert_eq!(app.active_help_copyable_text(), source);
    app.help_gui_markdown = "# Replaced\n\nNEW_BODY\n".into();
    let replaced = app.rendered_help_markdown_for(HelpDoc::Gui);
    assert!(!Arc::ptr_eq(&index, &app.ensure_gui_manual_pages()));
    assert_eq!(app.help_gui_page, 0);
    assert!(replaced.contains("NEW_BODY"));
    assert!(!replaced.contains("FIRST_BODY"));
}

#[test]
fn gui_manual_pages_search_whole_manual_and_navigate_to_matching_section() {
    let mut app = GENtleApp::default();
    app.help_gui_markdown = "# Manual\n\n## One\n\nneedle one\n\n## Two\n\nneedle two\n".into();
    app.help_search_query = "needle".into();
    app.refresh_help_search_matches();
    assert_eq!(app.help_search_matches.len(), 2);
    assert_eq!(app.help_gui_page, 1);
    assert!(
        app.help_source_markdown(HelpDoc::Gui)
            .contains("needle one")
    );
    app.select_next_help_match();
    assert_eq!(app.help_gui_page, 2);
    assert!(
        app.help_source_markdown(HelpDoc::Gui)
            .contains("needle two")
    );
    app.select_prev_help_match();
    assert_eq!(app.help_gui_page, 1);
    let before = serde_json::to_value(app.engine.read().unwrap().state()).unwrap();
    click_help_link(&mut app, "Next section");
    assert_eq!(app.help_gui_page, 2);
    assert_eq!(
        serde_json::to_value(app.engine.read().unwrap().state()).unwrap(),
        before
    );
}

#[test]
fn gui_manual_pages_real_help_wheel_moves_content_without_rendering_other_sections() {
    let mut app = GENtleApp::default();
    app.help_gui_markdown = format!(
        "# Manual\n\n## Scroll\n\n{}\n## Hidden\n\nHIDDEN_PAGE_BODY\n",
        (0..150)
            .map(|i| format!("Scrolling row {i}\n\n"))
            .collect::<String>()
    );
    app.select_gui_manual_page(1);
    app.show_help_dialog = true;
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let mut time = 0.0;
    // Exercise the actual hosted Help window. A synthetic CentralPanel on an
    // unregistered layer can paint text, but cannot own wheel input in egui.
    let mut frame = |events| {
        time += 1.0 / 60.0;
        let mut output = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1200.0, 1000.0),
                )),
                time: Some(time),
                events,
                ..Default::default()
            },
            |_ui| app.render_help_dialog(&ctx),
        );
        output.textures_delta.clear();
        output
    };
    let visible_rows = |output: &egui::FullOutput| {
        output
            .shapes
            .iter()
            .filter_map(|shape| {
                let egui::epaint::Shape::Text(text) = &shape.shape else {
                    return None;
                };
                assert!(!text.galley.job.text.contains("HIDDEN_PAGE_BODY"));
                if !egui::Rect::from_min_size(text.pos, text.galley.size())
                    .intersects(shape.clip_rect)
                {
                    return None;
                }
                text.galley
                    .job
                    .text
                    .strip_prefix("Scrolling row ")
                    .and_then(|s| s.trim().parse::<usize>().ok())
            })
            .min()
            .unwrap()
    };
    for _ in 0..3 {
        frame(vec![]).drop_without_applying_deltas();
    }
    let before = frame(vec![]);
    let first = visible_rows(&before);
    let pointer = before
        .shapes
        .iter()
        .find_map(|shape| help_link_text_center(&shape.shape, &format!("Scrolling row {first}")))
        .unwrap();
    before.drop_without_applying_deltas();
    let mut last = first;
    for index in 0..20 {
        let mut events = vec![egui::Event::PointerMoved(pointer)];
        if index == 0 {
            events.push(egui::Event::MouseWheel {
                unit: egui::MouseWheelUnit::Point,
                delta: egui::vec2(0.0, -350.0),
                modifiers: egui::Modifiers::NONE,
                phase: egui::TouchPhase::Move,
            });
        }
        let output = frame(events);
        last = visible_rows(&output);
        output.drop_without_applying_deltas();
    }
    assert!(
        last > first,
        "wheel must advance the visible manual rows: {first} -> {last}"
    );
    assert_eq!(app.help_gui_page, 1);
}

// Hand-written temporary Markdown and egui pointer events reproduce the links
// in 03.09 without rasterizing its historical screenshots or touching user files.
fn help_link_test_app(root: &Path, newline: &str) -> GENtleApp {
    let directory = root.join("help with ' spaces");
    fs::create_dir_all(&directory).unwrap();
    let pages = [
        (
            "03-09_synthetic_biology_design_build_verify.md",
            "Source",
            "[03.07 Gibson Arrangements](./03-07_gibson_arrangements_gui.md)\n\n[03.08 Gibson Physical Rack][rack]\n\n[rack]: ./03-08_gibson_physical_rack_gui.md\n",
        ),
        (
            "03-07_gibson_arrangements_gui.md",
            "Arrangements",
            "[Return](./03-09_synthetic_biology_design_build_verify.md)\n",
        ),
        (
            "03-08_gibson_physical_rack_gui.md",
            "Physical Rack",
            "Target content\n",
        ),
    ];
    let mut app = GENtleApp::default();
    for (name, title, body) in pages {
        let path = directory.join(name);
        fs::write(&path, format!("# {title}\n\n{body}").replace('\n', newline)).unwrap();
        app.help_tutorial_entries.push(
            GENtleApp::help_tutorial_entry_from_catalog_entry(catalog_entry(&path, title)).unwrap(),
        );
    }
    app.open_help_tutorial_doc(0);
    app
}

fn help_link_frame(
    app: &mut GENtleApp,
    ctx: &egui::Context,
    events: Vec<egui::Event>,
) -> egui::FullOutput {
    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(egui::Rect::from_min_size(
            egui::Pos2::ZERO,
            egui::vec2(1200.0, 1000.0),
        )),
        events,
        ..Default::default()
    });
    crate::egui_compat::show_central_panel_for_test_context(
        ctx,
        egui::CentralPanel::default(),
        |ui| app.render_help_contents(ui),
    );
    crate::egui_compat::end_test_pass(ctx)
}

fn help_link_text_center(shape: &egui::epaint::Shape, label: &str) -> Option<egui::Pos2> {
    match shape {
        egui::epaint::Shape::Text(text) if text.galley.job.text == label => {
            Some(text.pos + text.galley.size() * 0.5)
        }
        egui::epaint::Shape::Vec(shapes) => shapes
            .iter()
            .find_map(|shape| help_link_text_center(shape, label)),
        _ => None,
    }
}

fn click_help_link(app: &mut GENtleApp, label: &str) {
    let ctx = egui::Context::default();
    help_link_frame(app, &ctx, vec![]).drop_without_applying_deltas();
    let output = help_link_frame(app, &ctx, vec![]);
    let pos = output
        .shapes
        .iter()
        .find_map(|shape| help_link_text_center(&shape.shape, label))
        .expect("rendered Markdown hyperlink");
    output.drop_without_applying_deltas();
    for pressed in [true, false] {
        let output = help_link_frame(
            app,
            &ctx,
            vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ],
        );
        assert!(
            !output
                .platform_output
                .commands
                .iter()
                .any(|command| matches!(command, egui::OutputCommand::OpenUrl(_))),
            "a local chapter click must not launch a browser"
        );
        output.drop_without_applying_deltas();
    }
}

#[test]
fn help_markdown_links_click_inline_and_reference_chapters_without_project_changes() {
    let source =
        include_str!("../../../docs/tutorial/03-09_synthetic_biology_design_build_verify.md");
    for name in [
        "03-07_gibson_arrangements_gui.md",
        "03-08_gibson_physical_rack_gui.md",
    ] {
        assert!(source.contains(&format!("](./{name})")));
    }
    for newline in ["\n", "\r\n"] {
        let temp = tempfile::tempdir().unwrap();
        let mut app = help_link_test_app(temp.path(), newline);
        let before = serde_json::to_value(app.engine.read().unwrap().state()).unwrap();
        let revision = app.engine.read().unwrap().mutation_revision();
        click_help_link(&mut app, "03.07 Gibson Arrangements");
        assert_eq!(app.help_tutorial_selected, 1);
        assert_eq!(app.help_tutorial_title, "Arrangements");
        click_help_link(&mut app, "Return");
        assert_eq!(app.help_tutorial_selected, 0);
        click_help_link(&mut app, "03.08 Gibson Physical Rack");
        assert_eq!(app.help_tutorial_selected, 2);
        assert!(app.help_tutorial_markdown.contains("Target content"));
        assert_eq!(app.help_tutorial_entries.len(), 3, "reuse catalog entries");
        assert_eq!(app.help_doc, HelpDoc::Tutorial);
        assert!(app.show_help_dialog);
        assert!(app.help_link_notice.is_none());
        assert_eq!(
            serde_json::to_value(app.engine.read().unwrap().state()).unwrap(),
            before
        );
        assert_eq!(app.engine.read().unwrap().mutation_revision(), revision);
    }
}

#[test]
fn help_markdown_links_preserve_web_urls_clipboard_and_other_window_commands() {
    let temp = tempfile::tempdir().unwrap();
    let mut app = help_link_test_app(temp.path(), "\n");
    let ctx = egui::Context::default();
    ctx.begin_pass(egui::RawInput::default());
    ctx.open_url(egui::OpenUrl::same_tab("./other-window.md"));
    let start = ctx.output(|output| output.commands.len());
    ctx.open_url(egui::OpenUrl::same_tab(
        "./03-07_gibson_arrangements_gui.md",
    ));
    ctx.open_url(egui::OpenUrl::new_tab("https://example.invalid/guide.md"));
    ctx.copy_text("Keep clipboard output".into());
    app.follow_help_markdown_links(&ctx, start);
    assert_eq!(app.help_tutorial_selected, 1);
    ctx.output(|output| {
        assert_eq!(output.commands.len(), 3);
        assert!(matches!(&output.commands[0], egui::OutputCommand::OpenUrl(url) if url.url == "./other-window.md"));
        assert!(matches!(&output.commands[1], egui::OutputCommand::OpenUrl(url) if url.url == "https://example.invalid/guide.md" && url.new_tab));
        assert!(matches!(&output.commands[2], egui::OutputCommand::CopyText(text) if text == "Keep clipboard output"));
    });
    crate::egui_compat::discard_test_pass_output(&ctx);
}

#[test]
fn help_markdown_links_report_missing_targets_and_disclose_fragment_limit() {
    let temp = tempfile::tempdir().unwrap();
    let mut app = help_link_test_app(temp.path(), "\n");
    let ctx = egui::Context::default();
    ctx.begin_pass(egui::RawInput::default());
    ctx.open_url(egui::OpenUrl::same_tab("./missing.md"));
    app.follow_help_markdown_links(&ctx, 0);
    assert_eq!(app.help_tutorial_selected, 0);
    assert!(
        app.help_link_notice
            .as_deref()
            .unwrap()
            .contains("could not be found")
    );
    assert!(ctx.output(|output| output.commands.is_empty()));
    crate::egui_compat::discard_test_pass_output(&ctx);
    let notice = app.help_link_notice.clone().unwrap();
    let output = help_link_frame(&mut app, &ctx, vec![]);
    assert!(
        output
            .shapes
            .iter()
            .any(|shape| help_link_text_center(&shape.shape, &notice).is_some())
    );
    output.drop_without_applying_deltas();
    assert!(
        app.open_help_markdown_link("./03-07_gibson_arrangements_gui.md#details")
            .is_ok()
    );
    assert_eq!(app.help_tutorial_selected, 1);
    assert!(
        app.help_link_notice
            .as_deref()
            .unwrap()
            .contains("automatic section scrolling")
    );
    assert!(
        app.open_help_markdown_link("./03-08_gibson_physical_rack_gui.md")
            .is_ok()
    );
    assert!(app.help_link_notice.is_none());
}

#[test]
fn help_markdown_links_reuse_manual_tabs_and_discover_unlisted_documents() {
    let mut app = GENtleApp::default();
    let source = Path::new("docs/tutorial/03-09_synthetic_biology_design_build_verify.md");
    app.help_tutorial_entries.push(
        GENtleApp::help_tutorial_entry_from_catalog_entry(catalog_entry(source, "Source")).unwrap(),
    );
    app.help_doc = HelpDoc::Tutorial;
    app.open_help_markdown_link("../cli.md").unwrap();
    assert_eq!(app.help_doc, HelpDoc::Cli);
    assert!(app.show_help_dialog);
    let temp = tempfile::tempdir().unwrap();
    let mut app = help_link_test_app(temp.path(), "\n");
    let source_dir = app
        .help_document_path()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf();
    fs::write(source_dir.join("extra.md"), "# Extra linked document\n").unwrap();
    app.open_help_markdown_link("./extra.md").unwrap();
    assert_eq!(app.help_tutorial_title, "Extra linked document");
    assert_eq!(app.help_tutorial_entries.len(), 4);
    app.open_help_markdown_link("./extra.md").unwrap();
    assert_eq!(app.help_tutorial_entries.len(), 4);
}

#[test]
fn help_markdown_links_do_not_intercept_external_or_absolute_targets() {
    for destination in [
        "https://example.invalid/guide.md",
        "mailto:example@example.invalid",
        "file:///tmp/guide.md",
        "//example.invalid/guide.md",
        "/tmp/guide.md",
        "C:\\docs\\guide.md",
        "C:/docs/guide.md",
        "#section",
        "figure.svg",
    ] {
        assert!(
            GENtleApp::help_relative_markdown_link(destination).is_none(),
            "{destination}"
        );
    }
    assert_eq!(
        GENtleApp::help_relative_markdown_link("../guide.md#section"),
        Some("../guide.md")
    );
}

fn catalog_entry(path: &Path, title: &str) -> crate::workflow_examples::TutorialCatalogEntry {
    serde_json::from_value(serde_json::json!({
        "id": "synthetic-help-test",
        "title": title,
        "path": path.to_str().unwrap(),
        "type": "operational_reference",
        "status": "manual/reference",
        "source": "hand-written temporary test fixture",
        "audiences": ["agent_users"],
        "prerequisites": ["before-help-test"],
        "next_steps": ["after-help-test"],
        "group_label": "Getting Started & Interfaces",
        "group_order": 1,
        "group_position": 2,
        "decimal_id": "01.02"
    }))
    .unwrap()
}

fn assert_no_image_work(report: &serde_json::Value) {
    assert_eq!(report["dropped_events"], 0);
    assert_eq!(report["dropped_help_image_observations"], 0);
    assert_eq!(report["help_image_work"]["saturated"], false);
    for key in [
        "svg_references",
        "cache_hits",
        "preparation_failures",
        "rasterization_attempts",
        "rasterization_completed",
        "rasterization_failures",
        "rasterization_us",
    ] {
        assert_eq!(report["help_image_work"][key], 0, "{key}");
    }
}

#[test]
fn lazy_tutorial_heading_preserves_precedence_and_fallbacks_without_images() {
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("01_fallback_title.md");
    let (_, trace) = startup_trace::capture(|| {
        for newline in ["\n", "\r\n"] {
            let markdown = [
                "![Before heading](missing-before.svg)",
                "#",
                "## Full document heading, not the catalog abbreviation",
                "![After heading](missing-after.svg)",
            ]
            .join(newline);
            fs::write(&path, &markdown).unwrap();
            let entry = GENtleApp::help_tutorial_entry_from_catalog_entry(catalog_entry(
                &path,
                "Catalog abbreviation",
            ))
            .unwrap();
            assert_eq!(
                entry.title,
                "Full document heading, not the catalog abbreviation"
            );
            assert_eq!(entry.decimal_id.as_deref(), Some("01.02"));
            assert_eq!(entry.audiences, ["agent_users"]);
            assert_eq!(entry.prerequisites, ["before-help-test"]);
            assert_eq!(entry.next_steps, ["after-help-test"]);
            assert!(entry.summary.contains("status: manual/reference"));
            assert_eq!(fs::read_to_string(&path).unwrap(), markdown);
        }
        for contents in [b"No heading\n".as_slice(), b"# Invalid UTF-8\n\xff"] {
            fs::write(&path, contents).unwrap();
            for (catalog_title, expected) in [
                ("Catalog fallback", "Catalog fallback"),
                ("  ", "fallback title"),
            ] {
                let entry = GENtleApp::help_tutorial_entry_from_catalog_entry(catalog_entry(
                    &path,
                    catalog_title,
                ))
                .unwrap();
                assert_eq!(entry.title, expected);
            }
        }
        fs::remove_file(&path).unwrap();
        assert!(GENtleApp::load_help_tutorial_heading(&path).is_err());
        assert!(
            GENtleApp::help_tutorial_entry_from_catalog_entry(catalog_entry(&path, "Missing"))
                .is_none()
        );
    });
    assert_no_image_work(&trace);
}

#[test]
fn lazy_tutorial_heading_keeps_relative_images_as_raw_text_without_preparation() {
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("heading_images.md");
    let svg = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"2\" height=\"2\"><rect width=\"2\" height=\"2\"/></svg>";
    fs::write(temp.path().join("badge.svg"), svg).unwrap();
    let title = "Example ![Badge](badge.svg) and ![Unavailable](missing.svg)";

    let (_, trace) = startup_trace::capture(|| {
        for newline in ["\n", "\r\n"] {
            let markdown = format!("## {title}{newline}{newline}Tutorial body{newline}");
            fs::write(&path, &markdown).unwrap();
            let entry = GENtleApp::help_tutorial_entry_from_catalog_entry(catalog_entry(
                &path,
                "Catalog fallback",
            ))
            .unwrap();
            // Discovery must not leak rewritten file/cache URLs into titles.
            assert_eq!(entry.title, title);
            assert_eq!(fs::read_to_string(&path).unwrap(), markdown);
            assert_eq!(
                fs::read_to_string(temp.path().join("badge.svg")).unwrap(),
                svg
            );
        }
    });
    assert_no_image_work(&trace);
}

#[test]
fn tutorial_navigation_resolves_catalog_ids_to_compact_targets() {
    let mut app = GENtleApp::default();
    app.help_tutorial_entries = GENtleApp::discover_help_tutorial_entries();
    let selected = app
        .help_tutorial_entries
        .iter()
        .position(|entry| entry.tutorial_id == "tss_collection_gui")
        .expect("TSS collection tutorial");
    app.help_tutorial_selected = selected;
    let entry = &app.help_tutorial_entries[selected];

    let prerequisites = app.help_tutorial_navigation_targets(&entry.prerequisites);
    let next_steps = app.help_tutorial_navigation_targets(&entry.next_steps);
    assert_eq!(prerequisites.len(), 1);
    assert_eq!(prerequisites[0].1, "08.13");
    assert!(prerequisites[0].2.contains("PWM/PSSM"));
    assert_eq!(next_steps.len(), 1);
    assert_eq!(next_steps[0].1, "08.16");
    assert!(next_steps[0].2.contains("Inspect TFBS"));
}

#[test]
fn tutorial_navigation_switch_keeps_help_open_and_only_prepares_selected_page() {
    let temp = tempfile::tempdir().unwrap();
    let mut app = GENtleApp::default();
    for id in ["before", "after", "unselected"] {
        let path = temp.path().join(format!("{id}.md"));
        fs::write(&path, format!("# {id}\n\n![Example]({id}-missing.svg)\n")).unwrap();
        let mut entry = catalog_entry(&path, id);
        entry.id = id.to_string();
        entry.prerequisites = if id == "after" {
            vec!["before".into()]
        } else {
            vec![]
        };
        entry.next_steps = if id == "before" {
            vec!["after".into()]
        } else {
            vec![]
        };
        app.help_tutorial_entries
            .push(GENtleApp::help_tutorial_entry_from_catalog_entry(entry).unwrap());
    }
    app.open_help_tutorial_doc(0);
    app.pending_focus_viewports.clear();
    for (field, expected_index) in [("next_steps", 1), ("prerequisites", 0)] {
        let (targets, resolution) = startup_trace::capture(|| {
            let current = &app.help_tutorial_entries[app.help_tutorial_selected];
            let ids = if field == "next_steps" {
                &current.next_steps
            } else {
                &current.prerequisites
            };
            app.help_tutorial_navigation_targets(ids)
        });
        assert_no_image_work(&resolution);
        assert_eq!(targets.len(), 1);
        assert_eq!(targets[0].0, expected_index);
        let (_, switched) = startup_trace::capture(|| {
            assert!(app.set_help_tutorial_selected(targets[0].0));
        });
        assert_eq!(app.help_tutorial_selected, expected_index);
        assert_eq!(app.help_doc, HelpDoc::Tutorial);
        assert!(app.show_help_dialog);
        assert!(app.pending_focus_viewports.is_empty());
        assert_eq!(switched["help_image_work"]["svg_references"], 1);
        assert_eq!(switched["help_image_work"]["rasterization_failures"], 1);
        let id = &app.help_tutorial_entries[expected_index].tutorial_id;
        assert!(app.help_tutorial_markdown.starts_with(&format!("# {id}\n")));
        assert!(
            app.help_tutorial_markdown
                .contains(&format!("{id}-missing.svg"))
        );
        assert!(!app.help_tutorial_markdown.contains("unselected"));
    }
    assert!(
        app.help_tutorial_navigation_targets(&["missing".into()])
            .is_empty()
    );
}

#[test]
fn lazy_tutorial_discovery_and_project_menu_keep_public_catalog_without_images() {
    let catalog = load_tutorial_catalog(Path::new(DEFAULT_TUTORIAL_CATALOG_PATH)).unwrap();
    let ((entries, guided, fallback), trace) = startup_trace::capture(|| {
        let entries = GENtleApp::discover_help_tutorial_entries();
        let guided = GENtleApp::tutorial_project_guided_walkthrough_entries();
        let mut fallback = Vec::new();
        GENtleApp::ensure_agent_interfaces_tutorial_entry(&mut fallback);
        GENtleApp::ensure_agent_interfaces_tutorial_entry(&mut fallback);
        (entries, guided, fallback)
    });
    assert_no_image_work(&trace);
    assert_eq!(entries.len(), catalog.entries.len());
    assert_eq!(fallback.len(), 1);
    assert_eq!(fallback[0].title, AGENT_INTERFACES_TUTORIAL_TITLE);
    let mut expected_order = Vec::new();
    let mut expected_guided_order = Vec::new();
    for entry in &catalog.entries {
        let path = GENtleApp::resolve_runtime_doc_path(&entry.path).unwrap();
        let markdown = fs::read_to_string(&path).unwrap();
        let expected_title = GENtleApp::markdown_first_heading(&markdown).unwrap_or_else(|| {
            if entry.title.trim().is_empty() {
                GENtleApp::markdown_title_from_path(&path)
            } else {
                entry.title.clone()
            }
        });
        let row = entries
            .iter()
            .find(|row| Path::new(&row.path) == path)
            .unwrap();
        assert_eq!(row.title, expected_title);
        assert_eq!(row.group_label, entry.group_label);
        assert_eq!(row.group_order, entry.group_order);
        assert_eq!(row.group_position, entry.group_position);
        expected_order.push(row.clone());
        if entry.status == "manual/reference" {
            expected_guided_order.push(row.clone());
        }
    }
    GENtleApp::sort_help_tutorial_entries_by_audience_group(&mut expected_order);
    GENtleApp::sort_help_tutorial_entries_by_audience_group(&mut expected_guided_order);
    let identities = |rows: &[HelpTutorialDocEntry]| {
        rows.iter()
            .map(|row| (&row.path, &row.title, &row.summary))
            .map(|(path, title, summary)| (path.clone(), title.clone(), summary.clone()))
            .collect::<Vec<_>>()
    };
    assert_eq!(identities(&entries), identities(&expected_order));
    assert_eq!(identities(&guided), identities(&expected_guided_order));
}

#[test]
fn lazy_tutorial_open_prepares_only_selected_images_once_and_preserves_fallback() {
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("selected.md");
    let other = temp.path().join("other.md");
    let svg = temp.path().join("selected.svg");
    fs::write(
        &svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="2" height="2"><rect width="2" height="2" fill="#123456"/></svg>"##,
    )
    .unwrap();
    fs::write(&path, "# Selected\n\n![Example](selected.svg)\n").unwrap();
    fs::write(&other, "# Other\n\n![Missing](missing.svg)\n").unwrap();
    let png = GENtleApp::help_svg_png_cache_path(&svg.canonicalize().unwrap()).unwrap();
    let missing_svg = temp.path().join("missing.svg");
    let missing_png = GENtleApp::help_svg_png_cache_path(&missing_svg).unwrap();
    let mut app = GENtleApp::default();
    let (_, discovery) = startup_trace::capture(|| {
        app.help_tutorial_entries = vec![
            GENtleApp::help_tutorial_entry_from_catalog_entry(catalog_entry(&other, "Other"))
                .unwrap(),
        ];
    });
    assert_no_image_work(&discovery);
    assert!(!png.exists());
    let (_, opened) = startup_trace::capture(|| {
        app.open_help_tutorial_path(path.to_str().unwrap(), "Fallback", "Synthetic guide")
            .unwrap();
    });
    assert_eq!(opened["help_image_work"]["svg_references"], 1);
    assert_eq!(opened["help_image_work"]["rasterization_attempts"], 1);
    assert_eq!(opened["help_image_work"]["rasterization_completed"], 1);
    assert_eq!(opened["help_image_work"]["rasterization_failures"], 0);
    assert!(png.is_file());
    assert!(!missing_png.exists());
    let expected = format!(
        "# Selected\n\n![Example]({})\n",
        GENtleApp::file_uri_from_path(&png)
    );
    assert_eq!(app.help_tutorial_markdown, expected);
    assert_eq!(app.help_tutorial_title, "Selected");
    assert!(app.show_help_dialog);
    let selected = app.help_tutorial_selected;
    let (_, switched) = startup_trace::capture(|| {
        assert!(app.set_help_tutorial_selected(0));
    });
    assert_eq!(switched["help_image_work"]["svg_references"], 1);
    assert_eq!(switched["help_image_work"]["rasterization_failures"], 1);
    assert_eq!(
        app.help_tutorial_markdown,
        format!(
            "# Other\n\n![Missing]({})\n",
            GENtleApp::file_uri_from_path(&missing_svg)
        )
    );
    let (_, reopened) = startup_trace::capture(|| {
        assert!(app.set_help_tutorial_selected(selected));
    });
    assert_eq!(reopened["help_image_work"]["cache_hits"], 1);
    assert_eq!(reopened["help_image_work"]["rasterization_attempts"], 0);
    assert_eq!(app.help_tutorial_markdown, expected);
    fs::remove_file(png).unwrap();
}

#[test]
fn lazy_tutorial_open_keeps_unreadable_document_errors() {
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("unreadable.md");
    let mut app = GENtleApp::default();
    for bytes in [None, Some(b"# Invalid\n\xff".as_slice())] {
        if let Some(bytes) = bytes {
            fs::write(&path, bytes).unwrap();
        }
        let (result, trace) = startup_trace::capture(|| {
            app.open_help_tutorial_path(path.to_str().unwrap(), "Fallback", "Synthetic guide")
        });
        assert!(result.unwrap_err().contains("could not be loaded"));
        assert!(!app.show_help_dialog);
        assert!(app.help_tutorial_entries.is_empty());
        assert_no_image_work(&trace);
    }
}
