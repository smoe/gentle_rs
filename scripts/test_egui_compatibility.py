"""Offline GUI dependency/API guards, not Rust or native-capture acceptance.

Inspect committed manifests/source only. No Cargo, binaries, screenshots,
provider calls or dependency updates are performed by these tests.
"""

from pathlib import Path
import re
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read_toml(path):
    return tomllib.loads((ROOT / path).read_text(encoding="utf-8"))


def family(version):
    match = re.match(r"(\d+)\.(\d+)", version)
    if not match:
        raise ValueError(f"expected explicit GUI version family, got {version!r}")
    return tuple(map(int, match.groups()))


class EguiCompatibilityTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
