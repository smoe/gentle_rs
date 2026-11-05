"""Explicit like-library gene-contrast regressions."""

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
SPEC = importlib.util.spec_from_file_location("compare_rna_read_end_gene_contrasts", ROOT / "scripts" / "compare_rna_read_end_gene_contrasts.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)
HEADER = tuple(__import__("rank_rna_read_end_genes").HEADER)


def projection(path: Path, sample: str, center: int) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        for index, value in enumerate((center - 1, center, center + 1), 1):
            writer.writerow([sample, "sha256:" + "a" * 64, f"r{index}", "utr3_reach_bp", value, "observed", "resolved", "tx", "GENE"])


class ContrastTests(unittest.TestCase):
    def test_primary_pair_is_ranked_and_enriched_pairs_stay_directional(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            centers = {"barcode13": 100, "barcode14": 120, "barcode15": 100, "barcode16": 150, "barcode17": 100, "barcode18": 90}
            paths = []
            for sample, center in centers.items():
                path = root / f"{sample}.tsv"
                projection(path, sample, center)
                paths.append(path)
            contrasts = root / "contrasts.tsv"
            contrasts.write_text(
                "\t".join(MODULE.CONTRAST_HEADER) + "\n"
                "regular\tbarcode16\tbarcode15\tregular_cdna\tprimary\n"
                "enrich_a\tbarcode14\tbarcode13\tenrichment_a\tdirectional_countercheck\n"
                "enrich_b\tbarcode18\tbarcode17\tenrichment_b\tdirectional_countercheck\n",
                encoding="utf-8",
            )
            output = root / "out"
            MODULE.run(paths, contrasts, "Sample 2", output, 3, 0.25, 10)
            with (output / "gene_contrasts.tsv").open(encoding="utf-8", newline="") as handle:
                row = next(row for row in csv.DictReader(handle, delimiter="\t") if row["metric"] == "utr3_reach_bp")
            self.assertEqual(row["primary_delta_median_bp"], "50")
            self.assertEqual(json.loads(row["directional_deltas_bp_json"]), {"enrich_a": 20, "enrich_b": -10})
            self.assertEqual(row["directional_concordant_count"], "1")
            report = json.loads((output / "report.json").read_text(encoding="utf-8"))
            self.assertIn("never pooled", report["directional_policy"])

    def test_reused_barcode_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = []
            for sample in ("b13", "b14", "b15"):
                path = root / f"{sample}.tsv"
                projection(path, sample, 100)
                paths.append(path)
            contrast = root / "contrasts.tsv"
            contrast.write_text(
                "\t".join(MODULE.CONTRAST_HEADER) + "\n"
                "primary\tb14\tb13\tregular\tprimary\n"
                "again\tb14\tb15\tenriched\tdirectional_countercheck\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "reused across contrasts"):
                MODULE.summarize(paths, contrast, "Sample 2", 1, 0.25)


if __name__ == "__main__":
    unittest.main()
