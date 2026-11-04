"""Synthetic comparator regressions; see compare_rna_read_end_metrics_README.md."""

import csv
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "compare_rna_read_end_metrics.py"
SPEC = importlib.util.spec_from_file_location("read_end_comparison", SCRIPT)
COMPARISON = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARISON)


class CompareRnaReadEndMetricsTests(unittest.TestCase):
    def minimal_input(self, path, *, newline="\n", zero=False):
        text = "sample_id\tmetric\tvalue_bp\tobservation_status\n"
        for sample in ("cohort_A", "cohort_B"):
            for value in (0 if zero else 10, 20):
                text += f"{sample}\tend_reach\t{value}\tobserved\n"
            for status in ("left_censored", "right_censored", "missing"):
                text += f"{sample}\tend_reach\t\t{status}\n"
        path.write_bytes(text.replace("\n", newline).encode())

    def invoke(self, source, output, *args):
        return subprocess.run([sys.executable, str(SCRIPT), str(source), str(output),
                               "--quantile-points", "11", *args],
                              capture_output=True, text=True, encoding="utf-8")

    def write_fixture(self, path: Path) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample_id", "metric", "value_bp", "observation_status"])
            for metric_index, metric in enumerate(("read_length_bp", "utr5_reach_bp", "utr3_reach_bp")):
                for sample_index in range(6):
                    for value in (100, 200, 300, 400):
                        writer.writerow([f"barcode{13 + sample_index}", metric, value + 10 * metric_index + sample_index, "observed"])
                    writer.writerow([f"barcode{13 + sample_index}", metric, "", "missing"])

    def test_six_samples_three_metrics_produce_all_fifteen_pairs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / "metrics.tsv", root / "out"
            self.write_fixture(source)
            subprocess.run([sys.executable, str(SCRIPT), str(source), str(output)], check=True)
            report = json.loads((output / "report.json").read_text(encoding="utf-8"))
            self.assertEqual(report["schema"], "gentle.rna_read_metric_comparison.v1")
            self.assertEqual(report["pair_count_per_metric"], 15)
            self.assertEqual(len(report["pairs"]), 45)
            with (output / "qq_pairs.tsv").open(encoding="utf-8") as handle:
                self.assertEqual(sum(1 for _ in handle) - 1, 15 * 3 * 101)
            summary = (output / "summary.tsv").read_text(encoding="utf-8")
            self.assertIn("barcode13\tread_length_bp\t4\t0\t0\t1\t", summary)
            for metric in ("read_length_bp", "utr5_reach_bp", "utr3_reach_bp"):
                self.assertIn("gray diagonal is identity", (output / f"qq_{metric}.svg").read_text(encoding="utf-8"))

    def test_missing_or_censored_value_cannot_masquerade_as_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "bad.tsv"
            source.write_text(
                "sample_id\tmetric\tvalue_bp\tobservation_status\n"
                "a\tread_length_bp\t100\tobserved\n"
                "a\tread_length_bp\t0\tmissing\n",
                encoding="utf-8",
            )
            completed = subprocess.run([sys.executable, str(SCRIPT), str(source), str(root / "out")], text=True, capture_output=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("never encoded as zero", completed.stderr)

    def test_observed_zero_is_not_missing_or_censored(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / "input.tsv", root / "out"
            self.minimal_input(source, zero=True)
            completed = self.invoke(source, output)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = json.loads((output / "report.json").read_bytes())
            for counts in report["counts"]:
                self.assertEqual([counts[key] for key in
                                  ("observed", "left_censored", "right_censored", "missing")], [2, 1, 1, 1])
            self.assertEqual(report["pairs"][0]["points"][0], [0.0, 0.0, 0.0])
            self.assertEqual(report["value_policy"]["positive_metrics"], [])

    def test_every_unobserved_status_still_refuses_a_numeric_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for status in ("left_censored", "right_censored", "missing"):
                with self.subTest(status=status):
                    source = root / "bad.tsv"
                    source.write_text("sample_id\tmetric\tvalue_bp\tobservation_status\n"
                                      f"a\tmetric\t0\t{status}\n", encoding="utf-8")
                    completed = self.invoke(source, root / "out")
                    self.assertNotEqual(completed.returncode, 0)
                    self.assertIn("never encoded as zero", completed.stderr)
                    self.assertFalse((root / "out").exists())

    def test_positive_requirement_is_explicit_and_not_inferred_from_metric_name(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.tsv"
            self.minimal_input(source, zero=True)
            refused = self.invoke(source, root / "zero", "--positive-metric", "end_reach")
            self.assertNotEqual(refused.returncode, 0)
            self.assertIn("positive observed value", refused.stderr)
            self.assertFalse((root / "zero").exists())
            unknown = self.invoke(source, root / "unknown", "--positive-metric", "unlisted_metric")
            self.assertNotEqual(unknown.returncode, 0)
            self.assertIn("absent from the input", unknown.stderr)
            self.assertFalse((root / "unknown").exists())
            self.minimal_input(source)
            passed = self.invoke(source, root / "positive", "--positive-metric", "end_reach")
            self.assertEqual(passed.returncode, 0, passed.stderr)
            report = json.loads((root / "positive/report.json").read_bytes())
            self.assertEqual(report["value_policy"]["positive_metrics"], ["end_reach"])

    def test_type_seven_quantiles_have_an_independent_arithmetic_control(self):
        for probability, expected in ((0, 0), (0.25, 7.5), (0.5, 20), (0.75, 37.5), (1, 60)):
            self.assertEqual(COMPARISON.quantile([0, 10, 30, 60], probability), expected)

    def test_multiline_and_ragged_tsv_fail_without_mutation_under_lf_and_crlf(self):
        header = "sample_id\tmetric\tvalue_bp\tobservation_status\n"
        rows = ('a\tmetric\t"1\n2"\tobserved\n',
                '"a\nb"\tmetric\t10\tobserved\n',
                'a\tmetric\t"12\tobserved\n',
                "a\tmetric\t12\tobserved\textra\n", "a\tmetric\t12\n")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for newline in ("\n", "\r\n"):
                for index, row in enumerate(rows):
                    with self.subTest(newline=repr(newline), row=index):
                        source = root / "bad.tsv"
                        raw = (header + row).replace("\n", newline).encode()
                        source.write_bytes(raw)
                        completed = self.invoke(source, root / "out")
                        self.assertNotEqual(completed.returncode, 0)
                        self.assertNotIn("Traceback", completed.stderr)
                        self.assertIn("line", completed.stderr)
                        self.assertFalse((root / "out").exists())
                        self.assertEqual(source.read_bytes(), raw)

    def test_invalid_scalar_values_fail_before_output_reservation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for value in ("-1", "1.5", "1_000", str(COMPARISON.MAX_VALUE_BP + 1)):
                with self.subTest(value=value):
                    source = root / "bad.tsv"
                    source.write_text("sample_id\tmetric\tvalue_bp\tobservation_status\n"
                                      f"a\tmetric\t{value}\tobserved\n", encoding="utf-8")
                    completed = self.invoke(source, root / "out")
                    self.assertNotEqual(completed.returncode, 0)
                    self.assertNotIn("Traceback", completed.stderr)
                    self.assertFalse((root / "out").exists())

    def test_existing_directory_and_missing_parent_are_not_modified(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.tsv"
            self.minimal_input(source)
            output = root / "existing"
            output.mkdir()
            refused = self.invoke(source, output)
            self.assertNotEqual(refused.returncode, 0)
            self.assertEqual(list(output.iterdir()), [])
            sentinel = output / "summary.tsv"
            sentinel.write_bytes(b"owned by someone else\n")
            refused = self.invoke(source, output)
            self.assertNotEqual(refused.returncode, 0)
            self.assertEqual(sentinel.read_bytes(), b"owned by someone else\n")
            refused = self.invoke(source, root / "missing/out")
            self.assertNotEqual(refused.returncode, 0)
            self.assertFalse((root / "missing").exists())
            refused = self.invoke(source, root / "existing/../out")
            self.assertNotEqual(refused.returncode, 0)
            self.assertFalse((root / "out").exists())

    def test_racing_output_creation_preserves_foreign_bytes_and_has_no_report(self):
        for filename in ("summary.tsv", "qq_pairs.tsv", "qq_end_reach.svg", "report.json"):
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                source, output = root / "input.tsv", root / "out"
                self.minimal_input(source)
                raw = source.read_bytes()
                original_link = os.link
                sentinel = b"other writer owns this file\n"

                def competing_link(partial, destination):
                    if Path(destination).name == filename:
                        Path(destination).write_bytes(sentinel)
                    return original_link(partial, destination)

                with patch.object(COMPARISON.os, "link", side_effect=competing_link):
                    with self.assertRaises(FileExistsError):
                        COMPARISON.run(source, output, 11)
                self.assertEqual((output / filename).read_bytes(), sentinel)
                if filename != "report.json":
                    self.assertFalse((output / "report.json").exists())
                self.assertEqual(source.read_bytes(), raw)
                self.assertFalse(list(output.glob(".*.partial")))

    def test_dangling_output_symlink_is_refused_without_following_it(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output, outside = root / "input.tsv", root / "out", root / "outside.tsv"
            self.minimal_input(source)
            probe = root / "symlink-probe"
            try:
                probe.symlink_to(outside)
            except OSError as error:
                self.skipTest(f"Host does not permit symlinks: {error}")
            probe.unlink()
            original_link = os.link

            def competing_link(partial, destination):
                Path(destination).symlink_to(outside)
                return original_link(partial, destination)

            with patch.object(COMPARISON.os, "link", side_effect=competing_link):
                with self.assertRaises(FileExistsError):
                    COMPARISON.run(source, output, 11)
            self.assertFalse(outside.exists())
            self.assertTrue((output / "summary.tsv").is_symlink())
            self.assertFalse((output / "report.json").exists())

    def test_existing_output_symlink_is_not_reused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output, outside = root / "input.tsv", root / "out", root / "outside"
            self.minimal_input(source)
            try:
                output.symlink_to(outside, target_is_directory=True)
            except OSError as error:
                self.skipTest(f"Host does not permit symlinks: {error}")
            completed = self.invoke(source, output)
            self.assertNotEqual(completed.returncode, 0)
            self.assertFalse(outside.exists())
            self.assertTrue(output.is_symlink())

    def test_staging_handles_close_before_publication_and_report_is_last(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / "input.tsv", root / "out"
            self.minimal_input(source)
            original_link, published = os.link, []

            def inspect_link(partial, destination):
                partial = Path(partial)
                moved = partial.with_suffix(".closed")
                partial.rename(moved)
                moved.rename(partial)
                published.append(Path(destination).name)
                return original_link(partial, destination)

            with patch.object(COMPARISON.os, "link", side_effect=inspect_link):
                COMPARISON.run(source, output, 11)
            self.assertEqual(published[-1], "report.json")
            self.assertEqual(len(published), 4)
            self.assertFalse(list(output.glob(".*.partial")))
            for path in output.iterdir():
                moved = path.with_suffix(path.suffix + ".closed")
                path.rename(moved)
                moved.unlink()
            output.rmdir()

    def test_write_failure_retains_partial_run_without_a_completion_report(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / "input.tsv", root / "out"
            self.minimal_input(source)
            raw = source.read_bytes()
            with patch.object(COMPARISON, "svg_for_metric", side_effect=OSError("synthetic write failure")):
                with self.assertRaisesRegex(OSError, "synthetic write failure"):
                    COMPARISON.run(source, output, 11)
            self.assertTrue((output / "summary.tsv").is_file())
            self.assertFalse((output / "report.json").exists())
            self.assertEqual(source.read_bytes(), raw)
            self.assertFalse(list(output.glob(".*.partial")))
            with self.assertRaises(FileExistsError):
                COMPARISON.run(source, output, 11)

    def test_unavailable_hard_links_fail_closed_without_overwrite_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / "input.tsv", root / "out"
            self.minimal_input(source)
            with patch.object(COMPARISON.os, "link", side_effect=OSError("hard links unavailable")):
                with self.assertRaisesRegex(OSError, "hard links unavailable"):
                    COMPARISON.run(source, output, 11)
            self.assertEqual(list(output.iterdir()), [])

    def test_lf_crlf_and_quoted_unicode_paths_preserve_semantics_and_raw_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "metrics with 'quote' and \u00e4"
            root.mkdir()
            reports, sources, summaries = [], [], []
            for index, newline in enumerate(("\n", "\r\n")):
                source, output = root / f"input {index}.tsv", root / f"output {index}"
                self.minimal_input(source, newline=newline, zero=True)
                raw = source.read_bytes()
                completed = self.invoke(source, output)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                report = json.loads((output / "report.json").read_bytes())
                self.assertEqual(report["source"]["sha256"], hashlib.sha256(raw).hexdigest())
                self.assertEqual(source.read_bytes(), raw)
                sources.append(report.pop("source"))
                reports.append(report)
                summaries.append((output / "summary.tsv").read_bytes())
                self.assertNotIn(b"\r\n", (output / "qq_pairs.tsv").read_bytes())
                self.assertNotIn(b"\r\n", (output / "qq_end_reach.svg").read_bytes())
            self.assertEqual(reports[0], reports[1])
            self.assertEqual(summaries[0], summaries[1])
            self.assertNotEqual(sources[0]["sha256"], sources[1]["sha256"])


if __name__ == "__main__":
    unittest.main()
