"""Offline GUI dependency/API guards, not Rust or native-capture acceptance.

Inspect committed manifests/source only. No Cargo, binaries, screenshots,
provider calls or dependency updates are performed by these tests.
"""

from pathlib import Path
import re
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[1]
LEGACY_ID_API = re.compile(r"\bId::new\s*\(|\bui\.id\s*\(\s*\)")


def read_toml(path):
    return tomllib.loads((ROOT / path).read_text(encoding="utf-8"))


def family(version):
    match = re.match(r"(\d+)\.(\d+)", version)
    if not match:
        raise ValueError(f"expected explicit GUI version family, got {version!r}")
    return tuple(map(int, match.groups()))


class EguiCompatibilityTests(unittest.TestCase):
    def test_no_deprecated_id_apis_in_gui_sources(self):
        for directory in ("src", "crates", "vendor/egui_commonmark/src",
                          "vendor/egui_commonmark_backend/src"):
            for path in sorted((ROOT / directory).rglob("*.rs")):
                with self.subTest(path=str(path.relative_to(ROOT))):
                    match = LEGACY_ID_API.search(path.read_text(encoding="utf-8"))
                    self.assertIsNone(match, "use Id::unique or Ui::scope_id, preserving the salt")

    def test_deprecated_id_api_guard_distinguishes_current_apis(self):
        for legacy in ('egui::Id::new("stable")', 'Id::new(("stable", 1))', 'ui.id().with("row")'):
            self.assertIsNotNone(LEGACY_ID_API.search(legacy), legacy)
        for current in ('egui::Id::unique("stable")', 'ui.scope_id().with("row")',
                        'texture.id()', 'protocol.id()'):
            self.assertIsNone(LEGACY_ID_API.search(current), current)

    def test_gui_and_vendored_markdown_share_one_version_family(self):
        root = read_toml("Cargo.toml")
        expected = family(root["workspace"]["dependencies"]["egui"])
        for name in ("eframe", "egui_extras"):
            self.assertEqual(family(root["dependencies"][name]["version"]), expected, name)
        for path in ("vendor/egui_commonmark/Cargo.toml",
                     "vendor/egui_commonmark_backend/Cargo.toml"):
            manifest = read_toml(path)
            for name in ("egui", "egui_extras"):
                self.assertEqual(family(manifest["dependencies"][name]["version"]), expected, path)
            if "eframe" in manifest.get("dev-dependencies", {}):
                self.assertEqual(family(manifest["dev-dependencies"]["eframe"]["version"]), expected)
        packages = read_toml("Cargo.lock")["package"]
        for name in ("egui", "eframe", "egui_extras", "ecolor", "emath", "epaint"):
            versions = [row["version"] for row in packages if row["name"] == name]
            self.assertEqual(len(versions), 1, (name, versions))
            self.assertEqual(family(versions[0]), expected, name)

    def test_capture_uses_viewport_bound_one_shot_callbacks(self):
        source = (ROOT / "src/agent_help.rs").read_text()
        app = (ROOT / "src/app/routine_and_agent_assistant_ui.rs").read_text()
        for removed in ("egui::UserData", "egui::Event::Screenshot", "collect_egui_capture_events"):
            self.assertNotIn(removed, source + app)
        request = source.split("pub fn request_egui_viewport_capture_for(", 1)[1].split(
            "pub fn render_agent_help_button(", 1)[0]
        for token in ("ctx.send_viewport_cmd_to(", "viewport_id,", "ScreenshotCallback::new",
                      "claim_egui_capture(token.request_id)", "window_title: token.window_title",
                      "request_repaint_of(viewport_id)", "request_repaint_of(egui::ViewportId::ROOT)"):
            self.assertIn(token, request)
        self.assertIn("for event in take_capture_events()", app)

    def test_callback_regressions_preserve_capture_consent_and_replay_guards(self):
        source = (ROOT / "src/agent_help.rs").read_text()
        app_tests = (ROOT / "src/app/tests.rs").read_text()
        self.assertNotIn("egui::Event::Screenshot", app_tests)
        self.assertIn("callbacks[0].clone().complete(image.clone())", source)
        self.assertIn("callbacks[0].clone().complete(image)", source)
        self.assertIn("screenshot_callback.clone().complete(image.clone())", app_tests)
        self.assertIn("screenshot_callback.complete(image)", app_tests)
        for test in ("agent_screenshot_request_and_auto_run_emit_no_capture",
                     "declining_agent_screenshot_request_captures_nothing",
                     "approved_agent_screenshot_is_one_shot_and_previews_before_transport",
                     "unsupported_agent_screenshot_provider_refuses_before_capture",
                     "timed_out_agent_screenshot_capture_discards_late_reply_id"):
            self.assertIn(f"fn {test}()", app_tests)

    def test_capture_keeps_logical_target_without_obsolete_transport_field(self):
        app = (ROOT / "src/app.rs").read_text()
        capture = app.split("struct AgentScreenshotCapture {", 1)[1].split("\n}", 1)[0]
        self.assertNotIn("capture_viewport_id", capture)
        self.assertIn("target_viewport_id: ViewportId", capture)
        approval = (ROOT / "src/app/routine_and_agent_assistant_ui.rs").read_text().split(
            "pub(super) fn approve_agent_screenshot_request(", 1)[1].split(
            "fn agent_screenshot_followup_prompt(", 1)[0]
        self.assertIn("let capture_viewport_id = if ctx.embed_viewports()", approval)
        self.assertIn("target.viewport_id", approval)
        tests = (ROOT / "src/app/tests.rs").read_text()
        self.assertIn("fn agent_screenshot_targets_native_or_embedded_viewport_without_stored_transport_id()", tests)
        self.assertNotIn("pending.capture_viewport_id", tests)
        self.assertIsNone(re.search(r"(?<![\w:])ViewportId::", tests),
                          "app tests use the imported egui module, not an unimported ViewportId")

    def test_native_capture_fixture_registers_the_child_before_approval(self):
        tests = (ROOT / "src/app/tests.rs").read_text()
        fixture = tests.split(
            "fn agent_screenshot_targets_native_or_embedded_viewport_without_stored_transport_id()", 1
        )[1].split("\n#[test]", 1)[0]
        setup = fixture.split("app.approve_agent_screenshot_request(&ctx);", 1)[0]
        native = re.search(r"if !embedded\s*\{(.*?)\n\s*\}", setup, re.S)
        self.assertIsNotNone(native, "native test setup must keep the child viewport alive")
        self.assertIn("ctx.show_viewport_deferred(", native.group(1))
        self.assertIn("target,", native.group(1))
        self.assertIn("assert!(matches!(viewport.class, egui::ViewportClass::Deferred))", fixture,
                      "ViewportClass has no Debug implementation for assert_eq!")
        self.assertIn("viewport.viewport_ui_cb.is_some()", fixture)

    def test_scroll_clients_use_public_offsets_and_defer_edge_requests(self):
        menu = (ROOT / "src/app/tutorial_menus.rs").read_text()
        for path in ("src/app/tutorial_menus.rs", "src/app/main_lineage_ui.rs",
                     "src/app/rack_workspace_ui.rs"):
            source = (ROOT / path).read_text()
            self.assertNotIn(".state.offset", source)
            self.assertIn(".state.clamped_offset()", source)
        for token in ("remove_temp::<f32>(edge_scroll_id)", "area.vertical_scroll_offset(offset)",
                      "data.insert_temp(edge_scroll_id, offset)",
                      "fn tutorial_menu_pending_edge_scroll_yields_to_new_wheel_or_drag_input()"):
            self.assertIn(token, menu)

    def test_synthetic_wheel_events_declare_their_source(self):
        for path in ("src/app/help_docs/tests.rs", "src/app/tutorial_menus.rs",
                     "src/main_area_dna/tests.rs",
                     "crates/gentle-benchmarks/benches/gui_operations.rs"):
            source = (ROOT / path).read_text()
            events = re.findall(r"egui::Event::MouseWheel\s*\{(.*?)\n\s*\}", source, re.S)
            self.assertTrue(events, path)
            for index, event in enumerate(events):
                with self.subTest(path=path, event=index):
                    self.assertIn("source: egui::MouseWheelSource::Unknown", event)

    def test_edge_hover_request_key_is_shared_with_the_test_frame(self):
        menu = (ROOT / "src/app/tutorial_menus.rs").read_text()
        for token in ("fn edge_scroll_request_id(ui: &egui::Ui) -> egui::Id",
                      "let edge_scroll_id = edge_scroll_request_id(ui);",
                      "let pending_edge_scroll_id = edge_scroll_request_id(ui);",
                      "let pending_id = hover.pending_edge_scroll_id;",
                      "fn tutorial_menu_pending_edge_scroll_is_applied_on_idle_frame()"):
            self.assertTrue(token in menu, f"missing shared-key regression contract: {token}")
        self.assertFalse('hover.area.id.with("edge_scroll")' in menu,
                         "the ScrollArea output ID is not the independent pending-request key")

    def test_test_windows_keep_named_explicit_ids(self):
        source = (ROOT / "src/main_area_dna/tests.rs").read_text()
        for name in ("pcr_designer_worker_poll", "locus_inspector_test",
                     "rna_read_table_height_probe"):
            with self.subTest(name=name):
                self.assertTrue(f'egui::Area::new(egui::Id::unique("{name}"))' in source,
                                f"{name} needs an explicit stable ID")
        implicit = re.search(r'egui::Area::new\(\s*"[^"\n]+"\.into\(\)', source)
        self.assertIsNone(implicit, "test windows must not use removed string-to-ID conversion")


if __name__ == "__main__":
    unittest.main()
