"""Synthetic regressions for gene-wise RNA read-end ranking."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location("rank_rna_read_end_genes", ROOT / "scripts" / "rank_rna_read_end_genes.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def write_projection(path: Path, sample: str, gene_values: dict[str, dict[str, list[tuple[str, str]]]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(MODULE.HEADER)
        counter = 0
        for gene, metrics in gene_values.items():
            for metric, values in metrics.items():
                for status, value in values:
                    counter += 1
                    writer.writerow([sample, "sha256:" + "a" * 64, f"r{counter}", metric, value, status, "resolved", "tx1", gene])


class RankingTests(unittest.TestCase):
    def test_focus_sample_ranks_3prime_shift_and_keeps_controls(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            common = {
                "read_length_bp": [("observed", "100"), ("observed", "110")],
                "utr5_reach_bp": [("observed", "20"), ("observed", "22")],
            }
            write_projection(root / "sample2.tsv", "sample2", {"TP53": {**common, "utr3_reach_bp": [("observed", "80"), ("observed", "90")]}})
            write_projection(root / "sample1.tsv", "sample1", {"TP53": {**common, "utr3_reach_bp": [("observed", "40"), ("observed", "50")]}})
            write_projection(root / "sample3.tsv", "sample3", {"TP53": {**common, "utr3_reach_bp": [("observed", "50"), ("observed", "60")]}})
            output = root / "out"
            MODULE.run([root / "sample1.tsv", root / "sample2.tsv", root / "sample3.tsv"], "sample2", output, 2, 0.25, 10)
            with (output / "gene_summary.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            utr3 = next(row for row in rows if row["metric"] == "utr3_reach_bp")
            self.assertEqual(utr3["comparison_sample_count"], "2")
            self.assertEqual(utr3["median_sample_delta_bp"], "35.0")
            self.assertIn("TP53", (output / "top_genes.svg").read_text(encoding="utf-8"))
            report = json.loads((output / "report.json").read_text(encoding="utf-8"))
            self.assertIn("does not establish TP73", report["nonclaims"][0])

    def test_censoring_excludes_ranking_without_turning_bound_into_value(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metric = {"utr3_reach_bp": [("observed", "80"), ("right_censored", "90")]}
            write_projection(root / "sample2.tsv", "sample2", {"GENE": metric})
            write_projection(root / "sample1.tsv", "sample1", {"GENE": {"utr3_reach_bp": [("observed", "40"), ("observed", "50")]}})
            rows, _samples = MODULE.summarize([root / "sample1.tsv", root / "sample2.tsv"], "sample2", 1, 0.25)
            utr3 = next(row for row in rows if row["gene_id"] == "GENE" and row["metric"] == "utr3_reach_bp")
            self.assertEqual(utr3["focus_right_censored_n"], 1)
            self.assertEqual(utr3["ranking_status"], "insufficient_or_censored")


if __name__ == "__main__":
    unittest.main()
