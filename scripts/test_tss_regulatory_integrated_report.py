#!/usr/bin/env python3
"""Focused deterministic tests for the TSS-local integrated report helpers."""

from __future__ import annotations

import importlib.util
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))


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
APPEND = load_module(
    "append_tss_similarity_to_locus_report",
    "append_tss_similarity_to_locus_report.py",
)
BUNDLE = ROOT.parent / "docs" / "examples" / "regulatory_region_comparison" / \
    "tp73_tss_local_integrated"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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

    def test_receipt_bound_window_slice_is_assembly_forward(self) -> None:
        window = {"start_0based": "100", "end_0based_exclusive": "110", "strand": "+"}
        self.assertEqual(
            PREPARE.assembly_forward_slice("AACCGGTTAA", window, 103, 106),
            "CCGG",
        )
        assembly = "AACCGGTTAA"
        transcript_oriented = PREPARE.reverse_complement(assembly)
        window["strand"] = "-"
        self.assertEqual(
            PREPARE.assembly_forward_slice(transcript_oriented, window, 103, 106),
            "CCGG",
        )
        with self.assertRaises(RuntimeError):
            PREPARE.assembly_forward_slice("AACCGGTTAA", window, 99, 106)


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

    def test_reverse_chain_decreases_in_query_without_false_break(self) -> None:
        reverse = [
            {"sstart": 10, "send": 20, "qstart": 90, "qend": 80, "bitscore": 20},
            {"sstart": 30, "send": 40, "qstart": 60, "qend": 50, "bitscore": 20},
        ]
        self.assertEqual([broken for _, broken in APPEND.ordered_blocks(reverse)],
                         [False, False])
        reordered = [reverse[0], {**reverse[1], "qstart": 100, "qend": 95}]
        self.assertEqual([broken for _, broken in APPEND.ordered_blocks(reordered)],
                         [False, True])

    def test_similarity_panel_is_inserted_before_interpretation_footer(self) -> None:
        base = (
            '<svg height="100" viewBox="0 0 1400 100" width="1400">\n'
            '<rect fill="#ffffff" height="100" width="1400" x="0" y="0"/>\n'
            '<text x="34" y="20">existing evidence</text>\n'
            '<text data-gentle-overlay-non-claims="true" x="34" y="70">\n'
            'Reporter interpretation boundaries\n</text>\n'
            '<text x="34" y="90">Evidence provenance</text>\n</svg>\n'
        )
        candidates = {
            "stretches": [{
                "stretch_id": "G_tss_stretch_1", "gene": "G",
                "start_1based": 100, "end_1based": 800,
                "tss_windows": [{"strand": "+", "chromosome": "1", "tss_1based": 600}],
            }],
            "regions": [],
        }
        output, height = APPEND.append_section(
            base, "G", candidates, {}, {}, {},
            {"minimum_bp": 40, "minimum_identity": 80.0, "maximum_evalue": 1e-5},
        )
        self.assertGreater(height, 100)
        self.assertLess(output.index("TSS-local Ensembl-feature promoterome similarity"),
                        output.index("Reporter interpretation boundaries"))
        self.assertIn('data-gentle-shifted-footer="true"', output)
        self.assertIn(f'height="{height}"', output)


class RetainedBundleTests(unittest.TestCase):
    def test_manifest_and_derivative_receipts(self) -> None:
        entries = {}
        for line in (BUNDLE / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split("  ", 1)
            entries[name] = digest
        self.assertGreaterEqual(len(entries), 20)
        for name, digest in entries.items():
            self.assertEqual(sha256(BUNDLE / name), digest, name)
        for gene in ("CD44", "TGFB1", "SERPINE1"):
            stem = f"{gene}_luciferase_planning_EnsemblReg_15TF_with_TSS_similarity"
            receipt = json.loads((BUNDLE / f"{stem}.receipt.json").read_text())
            for name, expected in receipt["outputs"].items():
                self.assertEqual(f"sha256:{sha256(BUNDLE / name)}", expected, name)

    def test_retained_comparison_uses_corrected_policy_and_limit_audit(self) -> None:
        comparison = json.loads((BUNDLE / "comparison.json").read_text())
        self.assertEqual(comparison["schema"], RENDER.comparison_tools.SCHEMA)
        self.assertEqual(comparison["counting_policy_id"],
                         RENDER.comparison_tools.COUNTING_POLICY)
        queries = comparison["tasks"]["blastn"]["queries"]
        self.assertEqual(len(queries), 15)
        for row in queries:
            self.assertTrue(row["other_gene_exclusion_verified"])
            self.assertEqual(row["counts_are_lower_bounds"],
                             row["target_cap_reached"] or row["hsp_cap_reached"])

    def test_tall_reports_place_similarity_before_original_footer(self) -> None:
        for gene in ("CD44", "TGFB1", "SERPINE1"):
            svg = (BUNDLE / f"{gene}_luciferase_planning_EnsemblReg_15TF_with_TSS_similarity.svg").read_text()
            self.assertLess(svg.index("TSS-local Ensembl-feature promoterome similarity"),
                            svg.index("Reporter interpretation boundaries"))
            self.assertIn('data-gentle-shifted-footer="true"', svg)


if __name__ == "__main__":
    unittest.main()
