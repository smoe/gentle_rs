#!/usr/bin/env python3
"""Tests for provenance-bound tutorial screenshot publication."""

from __future__ import annotations

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts import publish_tutorial_gui_screenshots as publisher


class TutorialScreenshotPublicationTests(unittest.TestCase):
    def test_publish_preserves_capture_bytes_and_emits_verifiable_manifest(self) -> None:
        with TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            evidence = root / "evidence"
            checkpoints = evidence / "chapter" / "checkpoints"
            checkpoints.mkdir(parents=True)
            (evidence / "acceptance-report.json").write_text(
                '{"status":"pass"}\n', encoding="utf-8"
            )
            raw = checkpoints / "action.raw.png"
            raw.write_bytes(b"synthetic png")
            snapshot = checkpoints / "action.snapshot.json"
            snapshot.write_text('{"generation":1}\n', encoding="utf-8")
            for role in ("orientation", "context"):
                (checkpoints / f"action.{role}.svg").write_text(
                    '<svg><image href="action.raw.png"/></svg>\n', encoding="utf-8"
                )
            record = {
                "schema": publisher.EVIDENCE_SCHEMA,
                "source_revision": "a" * 40,
                "gentle_binary_sha256": "b" * 64,
                "chapter_id": "chapter",
                "prose_step": 2,
                "step_id": "action",
                "requested_target": {"semantic_id": "control.action"},
                "capture": {
                    "captured_at_unix_ms": 1,
                    "raw": {"sha256": publisher.sha256_file(raw)},
                },
                "semantic_snapshot": {
                    "canonical_sha256": "c" * 64,
                    "retained_file_sha256": publisher.sha256_file(snapshot),
                },
                "derived_views": [
                    {
                        "role": role,
                        "path": str(checkpoints / filename),
                        "sha256": publisher.sha256_file(checkpoints / filename),
                    }
                    for role, filename in (
                        ("orientation", "action.orientation.svg"),
                        ("interaction_context", "action.context.svg"),
                    )
                ],
            }
            (checkpoints / "action.screenshot.json").write_text(
                json.dumps(record), encoding="utf-8"
            )
            selection = root / "selection.json"
            selection.write_text(
                json.dumps(
                    {
                        "schema": publisher.SELECTION_SCHEMA,
                        "capture_date": "2026-09-08",
                        "source_revision": "a" * 40,
                        "chapters": [
                            {"chapter_id": "chapter", "checkpoints": ["action"]}
                        ],
                    }
                ),
                encoding="utf-8",
            )
            output = root / "published"

            publisher.publish(evidence, selection, output)
            publisher.check(output)

            self.assertEqual(
                raw.read_bytes(), (output / "chapter" / "action.raw.png").read_bytes()
            )
            manifest = publisher.load_json(output / "publication-manifest.json")
            self.assertEqual(manifest["source_revision"], "a" * 40)
            self.assertEqual(len(manifest["checkpoints"]), 1)


if __name__ == "__main__":
    unittest.main()
