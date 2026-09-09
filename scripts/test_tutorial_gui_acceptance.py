#!/usr/bin/env python3
"""Deterministic unit tests for the external tutorial GUI acceptance runner."""

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts import tutorial_gui_acceptance as acceptance


class TutorialGuiAcceptanceTests(unittest.TestCase):
    def test_sequence_oracle_identity_checks_bases_topology_ends_and_features(self) -> None:
        import copy
        project = {"sequences": {"s": {"seq": {
            "seq": [65, 67, 71, 84], "topology": "Linear", "features": [],
        }, "overhang": {}, "gc_content": {"cached": 1}}}}
        expected = acceptance.sequence_content_identity(project, "s")
        changed = copy.deepcopy(project)
        changed["sequences"]["s"]["gc_content"] = {"cached": 99}
        self.assertEqual(expected, acceptance.sequence_content_identity(changed, "s"))
        for field, value in [("seq", [65]), ("topology", "Circular"), ("features", ["changed"])]:
            changed = copy.deepcopy(project)
            changed["sequences"]["s"]["seq"][field] = value
            self.assertNotEqual(expected, acceptance.sequence_content_identity(changed, "s"))
        changed = copy.deepcopy(project)
        changed["sequences"]["s"]["overhang"] = {"forward_5": "A"}
        self.assertNotEqual(expected, acceptance.sequence_content_identity(changed, "s"))
        with self.assertRaises(acceptance.AcceptanceFailure):
            acceptance.sequence_content_identity(project, "missing")

    def test_canonical_contract_hash_ignores_object_key_order(self) -> None:
        left = {"steps": [{"id": "a", "subject": {"b": "2", "a": "1"}}]}
        right = {"steps": [{"subject": {"a": "1", "b": "2"}, "id": "a"}]}
        self.assertEqual(
            acceptance.sha256_bytes(acceptance.canonical_json_bytes(left)),
            acceptance.sha256_bytes(acceptance.canonical_json_bytes(right)),
        )

    def test_json_path_supports_array_indices(self) -> None:
        value = {"pairs": [{"forward": {"end": 42}}]}
        self.assertEqual(acceptance.json_path(value, "pairs.0.forward.end"), 42)
        with self.assertRaises(KeyError):
            acceptance.json_path(value, "pairs.1.forward.end")

    def test_fixed_shell_command_quotes_runner_owned_paths(self) -> None:
        command = acceptance.fixed_shell_command(
            ["facts", "eval", "@/tmp/run with spaces/fact.json"]
        )
        self.assertEqual(
            command,
            "facts eval '@/tmp/run with spaces/fact.json'",
        )

    def test_report_assertions_cover_values_relations_and_non_empty(self) -> None:
        report = {
            "schema": "gentle.primer_design_report.v1",
            "report_id": "demo",
            "roi_start": 20,
            "roi_end": 40,
            "pairs": [{"forward_end": 10, "reverse_start": 50}],
        }
        verifier = {
            "schema": "gentle.primer_design_report.v1",
            "required_fields": ["report_id", "pairs.0.forward_end"],
            "assertions": [
                {"kind": "value", "path": "report_id", "equals": "demo"},
                {"kind": "value", "path": "pairs", "non_empty": True},
                {
                    "kind": "relation",
                    "left_path": "pairs.0.forward_end",
                    "op": "lte",
                    "right_path": "roi_start",
                },
            ],
        }
        acceptance.assert_report_contract(report, verifier)
        report["pairs"] = []
        with self.assertRaises(acceptance.AcceptanceFailure):
            acceptance.assert_report_contract(report, verifier)

    def test_profile_selection_requires_typed_acceptance(self) -> None:
        manifest = {
            "chapters": [
                {"id": "a", "gui_acceptance": {"profile": "smoke"}},
                {"id": "b", "gui_acceptance": {"profile": "offline-core"}},
            ]
        }
        selected = acceptance.selected_chapters(manifest, [], "smoke")
        self.assertEqual([chapter["id"] for chapter in selected], ["a"])
        with self.assertRaises(acceptance.AcceptanceFailure):
            acceptance.selected_chapters(manifest, ["missing"], None)

    def test_view_only_contract_starts_from_a_satisfied_navigation_invariant(self) -> None:
        self.assertEqual(
            acceptance.expected_starter_completion_truth({"view_only": True}),
            "satisfied",
        )
        self.assertEqual(
            acceptance.expected_starter_completion_truth({}),
            "unsatisfied",
        )

    def test_environment_values_are_hashed_not_exposed(self) -> None:
        rows, names = acceptance.redacted_environment(
            {
                "PATH": "/usr/bin",
                "GENTLE_CACHE": "/private/cache",
                "MISTRAL_API_KEY": "secret-value",
            }
        )
        self.assertEqual(names, {"GENTLE_CACHE", "MISTRAL_API_KEY"})
        serialized = str(rows)
        self.assertNotIn("/private/cache", serialized)
        self.assertNotIn("secret-value", serialized)
        self.assertTrue(all(len(row["value_sha256"]) == 64 for row in rows))

    def test_pending_metadata_is_flushed_before_scientific_effect(self) -> None:
        self.assertTrue(
            acceptance.should_flush_pending_save_before_step(
                True, {"scientific_effect": True}
            )
        )
        self.assertFalse(
            acceptance.should_flush_pending_save_before_step(
                False, {"scientific_effect": True}
            )
        )
        self.assertFalse(
            acceptance.should_flush_pending_save_before_step(
                True, {"scientific_effect": False}
            )
        )

    def test_scoped_lookup_does_not_fall_back_unless_explicit(self) -> None:
        runner = object.__new__(acceptance.TutorialAcceptanceRun)
        snapshot = {
            "items": [
                {
                    "semantic_id": "window.dna_viewer",
                    "window_id": "window.dna_viewer",
                    "subject_scope": None,
                }
            ]
        }
        self.assertIsNone(
            runner.item_for(
                snapshot,
                "window.dna_viewer",
                subject_scope="subject-deadbeef",
            )
        )
        self.assertIsNotNone(
            runner.item_for(
                snapshot,
                "window.dna_viewer",
                subject_scope="subject-deadbeef",
                allow_unscoped_fallback=True,
            )
        )

    def test_screenshot_geometry_uses_physical_pixels_and_clamped_context(self) -> None:
        item = {
            "rect_logical_points": {
                "min_x": 100.0,
                "min_y": 50.0,
                "max_x": 180.0,
                "max_y": 90.0,
            },
            "pixels_per_point": 1.5,
        }
        rectangle = acceptance.semantic_pixel_rect(item)
        self.assertEqual(
            rectangle,
            {"min_x": 150, "min_y": 75, "max_x": 270, "max_y": 135},
        )
        crop = acceptance.padded_crop_rect(rectangle, 320, 200)
        self.assertEqual(crop["min_x"], 0)
        self.assertEqual(crop["min_y"], 0)
        self.assertEqual(crop["max_x"], 320)
        self.assertEqual(crop["max_y"], 200)

    def test_xwininfo_geometry_and_client_to_root_translation(self) -> None:
        native = acceptance.parse_xwininfo_geometry(
            """
  Absolute upper-left X:  1
  Absolute upper-left Y:  20
  Width: 800
  Height: 600
""",
            4194306,
        )
        self.assertEqual(
            native,
            acceptance.NativeWindowGeometry(4194306, 1, 20, 800, 600),
        )
        semantic_window = {
            "rect_logical_points": {
                "min_x": 0.0,
                "min_y": 0.0,
                "max_x": 800.0,
                "max_y": 600.0,
            },
            "pixels_per_point": 1.0,
        }
        target = {
            "rect_logical_points": {
                "min_x": 118.0,
                "min_y": 346.0,
                "max_x": 156.0,
                "max_y": 384.0,
            },
            "pixels_per_point": 1.0,
        }
        self.assertEqual(
            acceptance.native_screen_rect(target, semantic_window, native),
            {"min_x": 119, "min_y": 366, "max_x": 157, "max_y": 404},
        )

    def test_native_translation_rejects_mismatched_scale_and_bounds(self) -> None:
        native = acceptance.NativeWindowGeometry(7, 100, 200, 400, 300)
        semantic_window = {
            "rect_logical_points": {
                "min_x": 10.0,
                "min_y": 20.0,
                "max_x": 410.0,
                "max_y": 320.0,
            },
            "pixels_per_point": 1.0,
        }
        target = {
            "rect_logical_points": {
                "min_x": 20.0,
                "min_y": 30.0,
                "max_x": 40.0,
                "max_y": 50.0,
            },
            "pixels_per_point": 2.0,
        }
        with self.assertRaises(acceptance.AcceptanceFailure):
            acceptance.native_screen_rect(target, semantic_window, native)
        target["pixels_per_point"] = 1.0
        target["rect_logical_points"]["max_x"] = 500.0
        with self.assertRaises(acceptance.AcceptanceFailure):
            acceptance.native_screen_rect(target, semantic_window, native)

    def test_parent_network_namespace_is_explicit_and_distinct(self) -> None:
        acceptance.validate_parent_network_namespace("net:[22]", "net:[11]")
        for own, parent in [
            ("net:[11]", "net:[11]"),
            ("not-a-namespace", "net:[11]"),
            ("net:[22]", ""),
        ]:
            with self.assertRaises(acceptance.AcceptanceFailure):
                acceptance.validate_parent_network_namespace(own, parent)

    def test_ewmh_client_and_process_identifiers_are_parsed_exactly(self) -> None:
        self.assertEqual(
            acceptance.parse_ewmh_client_ids(
                "_NET_CLIENT_LIST(WINDOW): window id # 0x400002, 0x4001cb"
            ),
            [4194306, 4194763],
        )
        self.assertEqual(
            acceptance.parse_ewmh_process_id(
                "_NET_WM_PID(CARDINAL) = 3785622", 4194306
            ),
            3785622,
        )
        self.assertEqual(
            acceptance.parse_ewmh_client_ids(
                "_NET_CLIENT_LIST:  no such atom on any window."
            ),
            [],
        )
        self.assertIsNone(
            acceptance.parse_ewmh_process_id(
                "_NET_WM_PID:  not found.", 4194306
            )
        )
        with self.assertRaises(acceptance.AcceptanceFailure):
            acceptance.parse_ewmh_client_ids("unexpected")

    def test_png_dimensions_and_teaching_svg_bind_one_raw_capture(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            png = root / "step.raw.png"
            png.write_bytes(
                b"\x89PNG\r\n\x1a\n"
                + b"\x00\x00\x00\rIHDR"
                + (640).to_bytes(4, "big")
                + (480).to_bytes(4, "big")
            )
            self.assertEqual(acceptance.png_dimensions(png), (640, 480))
            output = root / "step.context.svg"
            acceptance.write_screenshot_view_svg(
                output,
                png,
                640,
                480,
                {"min_x": 100, "min_y": 80, "max_x": 500, "max_y": 380},
                {"min_x": 220, "min_y": 160, "max_x": 280, "max_y": 200},
                "chapter / step",
            )
            svg = output.read_text(encoding="utf-8")
            self.assertIn('href="step.raw.png"', svg)
            self.assertIn('x="-100"', svg)
            self.assertIn('width="60" height="40"', svg)
            self.assertIn("chapter / step", svg)


if __name__ == "__main__":
    unittest.main()
