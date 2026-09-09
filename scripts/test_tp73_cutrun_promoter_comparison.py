#!/usr/bin/env python3
"""Integrity checks for the retained TP73 promoterome comparison example."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import unittest


ROOT = (
    Path(__file__).resolve().parents[1]
    / "docs/examples/regulatory_region_comparison/tp73_cutrun_supported_first_comparison"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class Tp73PromoterComparisonExampleTests(unittest.TestCase):
    def test_candidates_retain_seven_supported_tss_windows(self) -> None:
        payload = json.loads((ROOT / "candidate_regions.json").read_text())
        self.assertEqual(
            payload["schema"], "gentle.regulatory_region_comparison_sequences.v1"
        )
        regions = payload["regions"]
        self.assertEqual(len(regions), 7)
        self.assertEqual(
            {row["gene_query"] for row in regions}, {"CD44", "TGFB1", "SERPINE1"}
        )
        for row in regions:
            evidence = row["cutrun_support"]["evidence"]
            self.assertTrue(
                any(
                    max(item["experimental_minus_matched_gfp"].values()) > 0
                    for item in evidence
                )
            )

    def test_summary_retains_complete_recurrence_accounting(self) -> None:
        with (ROOT / "summary.tsv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 7)
        self.assertTrue(all(row["target_cap_reached"] == "False" for row in rows))
        distal = next(
            row for row in rows
            if row["gene"] == "SERPINE1" and row["tss_1based"] == "101122158"
        )
        self.assertEqual(distal["other_genes_25pct"], "27563")
        self.assertEqual(distal["other_genes_50pct"], "0")

    def test_receipt_binds_retained_outputs(self) -> None:
        receipt = json.loads((ROOT / "receipt.json").read_text())
        self.assertFalse(receipt["target_cap_reached"])
        for name in ("candidate_regions.json", "candidate_regions.fa"):
            self.assertEqual(sha256(ROOT / name), receipt["inputs"][name], name)
        for name, expected in receipt["outputs"].items():
            self.assertEqual(sha256(ROOT / name), expected, name)


if __name__ == "__main__":
    unittest.main()
