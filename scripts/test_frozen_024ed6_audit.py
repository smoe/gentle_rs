"""Verify retained real audit records without GENtle, network or Rust execution.

Records originate in Glen's owner-delivered small bundle and the authorized CI
run at 024ed6c0. README.md distinguishes original records from log-derived JSON
and producer-reported hashes of omitted artifacts.
"""

import hashlib
import itertools
import json
from pathlib import Path
import statistics
import unittest


ROOT = Path(__file__).resolve().parents[1]
BASE = Path("docs/audits/glen_followups_20261008/candidate_024ed6c0")
SOURCE = "024ed6c0520d38c4221517219d44c792084de4e5"
LOCK = "ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b"


def retained_hashes():
    return {BASE / name: digest for digest, name in
            (line.split("  ", 1) for line in (ROOT / BASE / "SHA256SUMS").read_text().splitlines())}


def record(name):
    return json.loads((ROOT / BASE / name).read_bytes())


class FrozenAuditTests(unittest.TestCase):
    def test_original_record_hashes_and_subset_membership(self):
        hashes = retained_hashes()
        self.assertEqual(len(hashes), 14)
        for path, digest in hashes.items():
            with self.subTest(path=path):
                self.assertTrue(path.is_relative_to(BASE))
                self.assertNotIn("..", path.parts)
                self.assertEqual(hashlib.sha256((ROOT / path).read_bytes()).hexdigest(), digest)
        files = {path.relative_to(ROOT) for path in (ROOT / BASE).rglob("*") if path.is_file()}
        self.assertEqual(files, set(hashes) | {BASE / "README.md", BASE / "SHA256SUMS"})
        self.assertEqual(hashes[BASE / "glen/HANDOFF-SHA256SUMS"],
                         "4f147f98457b02186c883f18da308dc19b652ce38f724ca4fa4c1f4aa491fc90")

    def test_ci_receipts_preserve_source_and_unrun_acceptance(self):
        for name in ("native_gui.json", "scroll_receipt_from_log.json",
                     "agent_cli_receipt_from_log.json"):
            with self.subTest(name=name):
                receipt = record(name)
                self.assertEqual(receipt["source_revision"], SOURCE)
                self.assertEqual(receipt["cargo_lock_sha256"], LOCK)
                self.assertEqual(receipt["profile"], "dev")
                self.assertIs(receipt["package_accepted"], False)
        self.assertIs(record("agent_cli_receipt_from_log.json")["live_agent_accepted"], False)
        self.assertIs(record("native_gui.json")["human_scientific_approval"], False)

    def test_scroll_cases_are_complete_directional_and_content_bound(self):
        receipt = record("scroll_receipt_from_log.json")
        self.assertEqual(receipt["exit_code"], 0)
        self.assertIs(receipt["native_latency_measured"], False)
        observations = receipt["observations"]
        self.assertEqual(len(observations), 24)
        fixtures = {row["fixture"] for row in observations}
        screens = {row["screen"] for row in observations}
        self.assertEqual(len(fixtures), 3)
        self.assertEqual(len(screens), 4)
        self.assertEqual({(row["fixture"], row["screen"], row["delta_x_points"]) for row in observations},
                         set(itertools.product(fixtures, screens, (-6, 6))))
        content = {}
        for row in observations:
            before, after = row["before"], row["after"]
            self.assertEqual(before[1:], after[1:])
            self.assertTrue(0 <= after[0] <= after[2] - after[1])
            self.assertGreater((after[0] - before[0]) * -row["delta_x_points"], 0)
            self.assertIs(row["sequence_content_unchanged"], True)
            self.assertIs(row["native_latency_measured"], False)
            self.assertEqual(row["source_revision"], SOURCE)
            self.assertEqual(row["source_version"], "0.1.0-internal.12+git." + SOURCE)
            content.setdefault(row["fixture"], row["sequence_content_sha256"])
            self.assertEqual(content[row["fixture"]], row["sequence_content_sha256"])

    def test_native_guard_is_synthetic_and_matches_completion(self):
        ledger = record("native_ledger.json")
        self.assertEqual(ledger["status"], "pass")
        self.assertEqual(ledger["completion"]["truth"], "satisfied")
        self.assertEqual(ledger["preflight"]["starter_completion_truth"], "unsatisfied")
        self.assertEqual(ledger["preflight"]["oracle_completion_truth"], "satisfied")
        self.assertEqual(ledger["network_policy"], "offline")
        self.assertEqual(ledger["network_enforcement"], "linux_network_namespace")
        self.assertEqual(len(ledger["steps"]), 8)
        comparison = record("base_comparison.json")
        self.assertIs(comparison["synthetic"], True)
        for name in ("human_scientific_approval", "online_vkorc1_accepted", "native_screenshot"):
            self.assertIs(comparison[name], False)
        raw = comparison["raw_sequences"]
        self.assertEqual(raw["fragment"].upper(), raw["reference"].upper())
        differences = [(i, a, b) for i, (a, b) in
                       enumerate(zip(raw["reference"].upper(), raw["alternate"].upper())) if a != b]
        self.assertEqual(differences, [(6, "C", "T")])
        self.assertEqual({len(seq) for seq in raw.values()}, {comparison["length_bp"]})
        self.assertEqual(comparison["project_sha256"], ledger["final_project_sha256"])

    def test_comparison_statistics_and_distinct_source_classes(self):
        for comparison, expected_sources in (
            ("comparison-A", {"e723f46eb3fe58b15ba4bad5ff975024cd761437",
                              "64baba4f1cd21c9df731b232d9f82ce55c9d6495"}),
            ("comparison-B", {SOURCE}),
        ):
            with self.subTest(comparison=comparison):
                rows = record(comparison + "/runs.json")
                stats = record(comparison + "/statistics.json")
                self.assertEqual({row["source_sha"] for row in rows}, expected_sources)
                self.assertEqual(len(rows), stats["sample_count"])
                self.assertEqual(len(rows), 12)
                self.assertEqual(len({row["run"] for row in rows}), 12)
                self.assertFalse(any("fill" in row["run"] for row in rows))
                for group_name, group in stats["groups"].items():
                    selected = [row for row in rows if row["role"] + "/" + row["condition"] == group_name]
                    self.assertEqual(len(selected), group["n"])
                    self.assertEqual(group["n"], 3)
                    for key, metric in group["metrics"].items():
                        values = [float(row[key]) for row in selected if row.get(key) is not None]
                        if metric["status"] == "measured":
                            self.assertEqual(metric["values"], values)
                            self.assertEqual(metric["median"], statistics.median(values))
                            self.assertEqual(metric["min"], min(values))
                            self.assertEqual(metric["max"], max(values))
                self.assertTrue(all(row["status"] == "pass" and row["failed_spans"] == 0
                                    and row["dropped_events"] == 0 for row in rows))


if __name__ == "__main__":
    unittest.main()
