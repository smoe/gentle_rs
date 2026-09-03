#!/usr/bin/env python3
"""Deterministic unit tests for the external tutorial GUI acceptance runner."""

from __future__ import annotations

import unittest

from scripts import tutorial_gui_acceptance as acceptance


class TutorialGuiAcceptanceTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
