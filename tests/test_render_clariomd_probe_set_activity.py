"""Synthetic tests for the local Clariom D probe-set activity renderer.

All inputs are hand-crafted at test runtime: a three-row raw PM table, a tiny
SQLite feature map, and a three-row probe-set annotation table. No CEL-derived
or vendor data is used. The fixture builder below is the deterministic
recreation recipe and these tests are its only GENtle consumer.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import io
import json
import os
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "render_clariomd_probe_set_activity.py"
SAMPLES = [
    "P_SKMel29_AdGFP_1.CEL",
    "P_SKMel29_AdDNp73beta_1.CEL",
    "P_SKMel29_AdTAp73alpha_1.CEL",
    "P_SKMel29_AdGFP_2.CEL",
    "P_SKMel29_AdDNp73beta_2.CEL",
    "P_SKMel29_AdTAp73alpha_2.CEL",
    "P_SKMel29_AdGFP_3.CEL",
    "P_SKMel29_AdDNp73beta_3.CEL",
    "P_SKMel29_AdTAp73alpha_3.CEL",
]
PROBESET_COLUMNS = [
    "gene_symbols",
    "probeset_id",
    "seqname",
    "strand",
    "start",
    "stop",
    "probe_count",
    "transcript_cluster_id",
    "locus_type",
    "exon_id",
    "psr_id",
    "probeset_type",
    "psr_type",
    "junction_start_edge",
    "junction_stop_edge",
    "level",
    "has_cds",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_synthetic_inputs(root: Path) -> tuple[Path, Path, Path]:
    annotation = root / "probe_sets.tsv"
    annotation_rows = [
        ("TP73", "PS_TP73_LOW", "200", "101"),
        ("TP73", "PS_TP73_HIGH", "100", "102"),
        ("FUS", "PS_FUS", "300", "201"),
    ]
    with annotation.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=PROBESET_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for gene, probeset_id, start, _probe_id in annotation_rows:
            writer.writerow(
                {
                    "gene_symbols": gene,
                    "probeset_id": probeset_id,
                    "seqname": "chr1" if gene == "TP73" else "chr16",
                    "strand": "+",
                    "start": start,
                    "stop": str(int(start) + 24),
                    "probe_count": "1",
                    "transcript_cluster_id": f"TC_{gene}",
                    "locus_type": "single",
                    "exon_id": f"EX_{probeset_id}",
                    "psr_id": f"PSR_{probeset_id}",
                    "probeset_type": "main->psr",
                    "psr_type": "main",
                    "level": "core",
                    "has_cds": "1",
                }
            )

    database = root / "platform.sqlite"
    connection = sqlite3.connect(database)
    try:
        connection.executescript(
            """
            CREATE TABLE featureSet (fsetid INTEGER PRIMARY KEY, man_fsetid TEXT NOT NULL);
            CREATE TABLE pmfeature (fsetid INTEGER NOT NULL, fid INTEGER NOT NULL, x INTEGER, y INTEGER);
            """
        )
        for fsetid, (_gene, probeset_id, _start, probe_id) in enumerate(
            annotation_rows, start=1
        ):
            connection.execute(
                "INSERT INTO featureSet(fsetid, man_fsetid) VALUES (?, ?)",
                (fsetid, probeset_id),
            )
            connection.execute(
                "INSERT INTO pmfeature(fsetid, fid, x, y) VALUES (?, ?, ?, ?)",
                (fsetid, int(probe_id), fsetid, fsetid + 10),
            )
        connection.commit()
    finally:
        connection.close()

    raw_features = root / "raw_features.tsv"
    values = {
        # log2(value + 1): GFP=2, DN=3, TA=4; paired contrasts are +1 and +2.
        101: {"GFP": 3, "DNp73beta": 7, "TAp73alpha": 15},
        # Paired contrasts are +3 and +4; the TP73 median must therefore be +2/+3.
        102: {"GFP": 3, "DNp73beta": 31, "TAp73alpha": 63},
        # Paired contrasts are -1 and 0.
        201: {"GFP": 7, "DNp73beta": 3, "TAp73alpha": 7},
    }
    with raw_features.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["probe_id", "x", "y", *SAMPLES],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for probe_id, condition_values in values.items():
            row: dict[str, int | str] = {"probe_id": probe_id, "x": 1, "y": 2}
            for sample in SAMPLES:
                condition = next(
                    name for name in ("DNp73beta", "GFP", "TAp73alpha") if name in sample
                )
                row[sample] = condition_values[condition]
            writer.writerow(row)
    return raw_features, database, annotation


class ClariomProbeSetActivityRendererTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.raw_features, self.database, self.annotation = write_synthetic_inputs(
            self.root
        )

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def run_renderer(self, output: Path, hash_seed: str) -> subprocess.CompletedProcess[str]:
        environment = os.environ.copy()
        environment["MPLCONFIGDIR"] = str(self.root / "matplotlib-config")
        environment["PYTHONHASHSEED"] = hash_seed
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--raw-features",
                str(self.raw_features),
                "--sqlite",
                str(self.database),
                "--fixture-probesets",
                str(self.annotation),
                "--vendor-probeset-zip",
                str(self.root / "unused-vendor.zip"),
                "--output-dir",
                str(output),
                "--genes",
                "TP73,FUS",
            ],
            cwd=ROOT,
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_repeat_renders_are_byte_identical_and_pairing_is_explicit(self) -> None:
        output_a = self.root / "render-a"
        output_b = self.root / "render-b"
        first = self.run_renderer(output_a, "101")
        second = self.run_renderer(output_b, "202")
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertEqual(second.returncode, 0, second.stderr)

        files_a = sorted(path.name for path in output_a.iterdir() if path.is_file())
        files_b = sorted(path.name for path in output_b.iterdir() if path.is_file())
        self.assertEqual(files_a, files_b)
        self.assertIn("gene_contrast_probe_set_summary.svg", files_a)
        self.assertIn("gene_contrast_probe_set_summary.pdf", files_a)
        self.assertIn("probe_set_individual_arrays_heatmap.pdf", files_a)
        self.assertNotIn("probe_set_individual_arrays_heatmap_10_gene.pdf", files_a)
        for name in files_a:
            self.assertEqual(sha256(output_a / name), sha256(output_b / name), name)

        manifest = json.loads((output_a / "manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["schema"], "gentle.clariomd_probe_set_activity.v1")
        self.assertEqual(
            [pairing["pair_id"] for pairing in manifest["sample_pairings"]],
            ["E1", "E2", "E3"],
        )
        self.assertEqual(
            manifest["sample_pairings"][0],
            {
                "pair_id": "E1",
                "GFP": "P_SKMel29_AdGFP_1.CEL",
                "DNp73beta": "P_SKMel29_AdDNp73beta_1.CEL",
                "TAp73alpha": "P_SKMel29_AdTAp73alpha_1.CEL",
            },
        )
        self.assertEqual(manifest["samples"], SAMPLES)
        self.assertIn(
            "matched GFP",
            manifest["contrast_definitions"]["paired_condition_minus_GFP"],
        )
        self.assertTrue(
            any("not significance tests" in item for item in manifest["interpretation_limitations"])
        )
        self.assertEqual(len(manifest["inputs"]), 3)
        self.assertTrue(all(len(item["sha256"]) == 64 for item in manifest["inputs"]))
        readme = (output_a / "README.md").read_text(encoding="utf-8")
        self.assertNotIn("unused-vendor.zip", readme)
        self.assertIn(manifest["inputs"][0]["sha256"], readme)

        with (output_a / "paired_gene_level_summary.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            paired_rows = list(csv.DictReader(handle, delimiter="\t"))
        tp73 = next(row for row in paired_rows if row["primary_gene"] == "TP73")
        self.assertEqual(tp73["E1_median_DN_minus_GFP"], "2.0000")
        self.assertEqual(tp73["E1_median_TA_minus_GFP"], "3.0000")

        with (output_a / "probe_set_activity_summary.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle, delimiter="\t"))
        low = next(row for row in summary_rows if row["probeset_id"] == "PS_TP73_LOW")
        self.assertEqual(low["log2_DNp73beta_minus_GFP"], "1.0000")
        self.assertEqual(low["log2_TAp73alpha_minus_GFP"], "2.0000")

        svg = (output_a / "probe_set_activity_heatmap.svg").read_text(encoding="utf-8")
        self.assertIn("Abundance: log2(raw PM mean + 1)", svg)
        self.assertIn("Contrasts: log2 difference", svg)
        self.assertIn("not a formal expression model", svg)

        plot_svg = (output_a / "gene_contrast_probe_set_summary.svg").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("<dc:date>", plot_svg)
        pdf = (output_a / "gene_contrast_probe_set_summary.pdf").read_bytes()
        self.assertIn(b"CreationDate", pdf)
        self.assertIn(b"1970", pdf)

    def test_bounded_line_reader_and_duplicate_gene_failure_are_useful(self) -> None:
        spec = importlib.util.spec_from_file_location("clariom_renderer", SCRIPT)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with self.assertRaisesRegex(RuntimeError, "line limit"):
            list(
                module.bounded_lines(
                    io.StringIO("header\nrow\n"),
                    source=Path("synthetic.tsv"),
                    max_lines=1,
                )
            )

        failed = subprocess.run(
            [sys.executable, str(SCRIPT), "--genes", "TP73,TP73"],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("--genes contains duplicate symbols: TP73", failed.stderr)

    def test_invalid_raw_intensity_names_the_sample_and_leaves_no_output(self) -> None:
        with self.raw_features.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = reader.fieldnames
            rows = list(reader)
        self.assertIsNotNone(fieldnames)
        rows[0]["P_SKMel29_AdGFP_1.CEL"] = "-1"
        with self.raw_features.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)

        output = self.root / "failed-render"
        failed = self.run_renderer(output, "303")
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("Raw intensity must be finite and non-negative", failed.stderr)
        self.assertIn("P_SKMel29_AdGFP_1.CEL", failed.stderr)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
