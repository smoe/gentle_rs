#!/usr/bin/env python3
"""Focused deterministic tests for the TSS-local integrated report helpers."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parent


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PREPARE = load_module(
    "prepare_tss_regulatory_similarity_candidates",
    "prepare_tss_regulatory_similarity_candidates.py",
)
RENDER = load_module(
    "render_integrated_tss_regulatory_report",
    "render_integrated_tss_regulatory_report.py",
)


class TssWindowTests(unittest.TestCase):
    def test_transcript_oriented_windows(self) -> None:
        self.assertEqual(PREPARE.selected_window(1_000, "+", 500, 200), (500, 1_200))
        self.assertEqual(PREPARE.selected_window(1_000, "-", 500, 200), (800, 1_500))
        with self.assertRaises(ValueError):
            PREPARE.selected_window(1_000, ".", 500, 200)

    def test_connected_stretches_merge_only_touching_windows(self) -> None:
        windows = [
            {"start_1based": 900, "end_1based": 1_200, "name": "b"},
            {"start_1based": 500, "end_1based": 900, "name": "a"},
            {"start_1based": 1_202, "end_1based": 1_500, "name": "c"},
        ]
        stretches = PREPARE.connected_stretches(windows)
        self.assertEqual([(row["start_1based"], row["end_1based"])
                          for row in stretches], [(500, 1_200), (1_202, 1_500)])
        self.assertEqual([row["name"] for row in stretches[0]["tss_windows"]], ["a", "b"])


class FrequencyTests(unittest.TestCase):
    def test_frequency_segments_count_distinct_genes_not_transcripts(self) -> None:
        rows = [
            {
                "promoter_id": "p1",
                "gene_ids": "g1",
                "query_intervals_0based_half_open": "0-10;20-30",
            },
            {
                "promoter_id": "p2",
                "gene_ids": "g1;g2",
                "query_intervals_0based_half_open": "5-25",
            },
        ]
        self.assertEqual(
            RENDER.frequency_segments(rows),
            [(0, 5, 1), (5, 10, 2), (10, 20, 2), (20, 25, 2), (25, 30, 1)],
        )

    def test_interval_merge_is_half_open(self) -> None:
        self.assertEqual(
            RENDER.merge_intervals([(8, 12), (0, 5), (5, 9), (20, 21)]),
            [(0, 12), (20, 21)],
        )


if __name__ == "__main__":
    unittest.main()
