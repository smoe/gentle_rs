//! Lazy Help regressions using the public catalog and temporary, hand-written
//! Markdown/SVG fixtures recreated below. No user files or font settings change.

use super::*;
use crate::gui_profiler::startup_trace;

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
