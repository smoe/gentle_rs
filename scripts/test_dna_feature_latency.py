"""Synthetic receipts/logs only; no Rust build, GUI, network or private fixtures.

Recreate with python3 -m unittest scripts.test_dna_feature_latency -v.
Exercises the prebuilt audit runner's identity, isolation and failure boundaries.
"""

import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from scripts import dna_feature_latency as audit


def observations(revision, workload=audit.LEGACY_WORKLOAD):
    fixtures = [(length, count) for length in (20_000, 250_000, 2_000_000)
                for count in (100, 1_000, 10_000)]
    if workload == audit.BOUNDARY_WORKLOAD:
        fixtures.append((250_000, 5_000))
    for length, count in fixtures:
        for interaction in ("steady", "pan_1bp", "zoom", "toggle_mrna", "select", "hover",
                            "resize_compact", "resize_desktop", "resize_fullhd"):
            row = {
                "schema": "gentle.dna_feature_latency_observation.v1", "revision": revision,
                "fixture_sha256": "f" * 64, "length_bp": length, "feature_count": count,
                "interaction": interaction,
                "before": {"tree_builds": 1, "layer_builds": 1, "layer_gc_bases": length, "linear": {"layouts": 1}},
                "after": {"tree_builds": 2, "layer_builds": 2, "layer_gc_bases": 2 * length, "linear": {"layouts": 2}},
            }
            if workload == audit.BOUNDARY_WORKLOAD:
                row["workload"] = workload
                row["derived_layers"] = (length, count) == (250_000, 5_000)
                for side in ("layers_before", "layers_after"):
                    row[side] = {
                        "annotation_features": count, "viewport": [0, 5_000],
                        "gc_bin_size_bp": 100,
                        **{name: {"total": 2, "viewport_eligible": 1, "enabled": True}
                           for name in ("restriction_groups", "gc_bins", "orfs", "methylation_sites")},
                    }
            yield row


class DnaFeatureLatencyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="DNA audit space ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.binary = self.root / "synthetic binary"
        self.binary.write_bytes(b"synthetic executable identity, never launched")
        self.receipt = self.root / "build.json"
        self.revision = "a" * 40
        audit.write_json(self.receipt, {
            "schema": "gentle.dna_feature_latency_build.v1", "profile": "dev",
            "binary": str(self.binary), "binary_sha256": audit.sha256(self.binary),
            "source": {"revision": self.revision, "dirty": False},
        })
        self.identity = {"schema": "gentle.dna_feature_latency_binary.v1", "revision": self.revision}

    def fake_run(self, argv, **kwargs):
        self.assertEqual(argv[0], str(self.binary))  # Never Cargo or a shell.
        self.assertNotIn("shell", kwargs)
        if argv[1:] == ["--gentle-identity"]:
            return subprocess.CompletedProcess(argv, 0, json.dumps(self.identity).encode(), b"")
        self.assertEqual(argv[1:], ["--test"])
        self.assertNotIn("GENTLE_PRIVATE_RESOURCE", kwargs["env"])
        self.assertEqual(kwargs["env"]["RUST_MIN_STACK"], "16777216")
        self.assertTrue(Path(kwargs["env"]["HOME"]).is_dir())
        self.assertTrue(str(kwargs["cwd"]) in kwargs["env"]["CRITERION_HOME"])
        for row in observations(self.revision, self.identity.get("workload", audit.LEGACY_WORKLOAD)):
            kwargs["stdout"].write((audit.PREFIX + json.dumps(row) + "\n").encode())
        return subprocess.CompletedProcess(argv, 0)

    def test_prebuilt_smoke_binds_logs_isolates_environment_and_never_builds(self):
        with patch.dict(audit.os.environ, {"GENTLE_PRIVATE_RESOURCE": "ignored"}), patch.object(audit.subprocess, "run", side_effect=self.fake_run) as run:
            result = audit.run(self.receipt, self.root / "run", "smoke", 10)
        self.assertEqual(run.call_count, 2)
        self.assertIsNone(result["failure"])
        self.assertFalse(result["native_gui_measured"])
        self.assertFalse(result["interactive_boundary_exercised"])
        self.assertEqual(len(result["observations"]), 81)
        self.assertEqual(result["log_sha256"], audit.sha256(self.root / "run/run.log"))
        self.assertEqual(result["runner_sha256"], audit.sha256(Path(audit.__file__)))
        self.assertEqual(len((self.root / "run/work.tsv").read_text().splitlines()), 82)

    def test_prepare_binds_cargo_artifact_and_rejects_source_change(self):
        def build(argv, **kwargs):
            if "--gentle-identity" in argv:
                return self.fake_run(argv, **kwargs)
            self.assertEqual(argv[:4], ["cargo", "test", "--locked", "--offline"])
            self.assertIn("--no-run", argv)
            self.assertNotIn("shell", kwargs)
            event = {"reason": "compiler-artifact", "target": {"name": "dna_feature_latency"},
                     "executable": str(self.binary)}
            kwargs["stdout"].write((json.dumps(event) + "\n").encode())
            return subprocess.CompletedProcess(argv, 0)

        source = {"revision": self.revision}
        with patch.object(audit.subprocess, "run", side_effect=build), patch.object(audit, "command_text", return_value="synthetic toolchain"), patch.object(audit, "source_identity", return_value=source):
            record = audit.prepare(self.root, self.root / "prepared", "dev", 10)
        self.assertEqual(record["binary_sha256"], audit.sha256(self.binary))
        self.assertEqual(record["source"], source)
        with patch.object(audit.subprocess, "run", side_effect=build), patch.object(audit, "source_identity", side_effect=[source, {"revision": "b" * 40}]), self.assertRaisesRegex(ValueError, "source changed"):
            audit.prepare(self.root, self.root / "changed", "dev", 10)
        self.assertFalse((self.root / "changed/build.json").exists())

    def test_audit_cannot_succeed_without_complete_statistics(self):
        build = json.loads(self.receipt.read_bytes())
        build["profile"] = "bench-audit"
        audit.write_json(self.receipt, build)

        def run(argv, **kwargs):
            if "--gentle-identity" in argv:
                return self.fake_run(argv, **kwargs)
            self.assertEqual(argv[1:], ["--bench", "--quick", "--noplot"])
            for row in observations(self.revision):
                kwargs["stdout"].write((audit.PREFIX + json.dumps(row) + "\n").encode())
            return subprocess.CompletedProcess(argv, 0)

        with patch.object(audit.subprocess, "run", side_effect=run), self.assertRaisesRegex(ValueError, "117 Criterion"):
            audit.run(self.receipt, self.root / "audit", "audit", 10)

    def test_tampered_binary_rejected_before_execution(self):
        self.binary.write_bytes(b"changed")
        with patch.object(audit.subprocess, "run") as run, self.assertRaisesRegex(ValueError, "binary changed"):
            audit.run(self.receipt, self.root / "run", "smoke", 10)
        run.assert_not_called()
        self.assertFalse((self.root / "run").exists())

    def test_wrong_compiled_revision_is_rejected(self):
        self.identity["revision"] = "b" * 40
        with patch.object(audit.subprocess, "run", side_effect=self.fake_run), self.assertRaisesRegex(ValueError, "identity"):
            audit.load_build(self.receipt, "smoke")

    def test_timing_rejects_development_build(self):
        with patch.object(audit.subprocess, "run") as run, self.assertRaisesRegex(ValueError, "bench-audit"):
            audit.load_build(self.receipt, "audit")
        run.assert_not_called()

    def test_timing_rejects_dirty_source(self):
        build = json.loads(self.receipt.read_bytes())
        build["profile"] = "bench-audit"
        build["source"]["dirty"] = True
        audit.write_json(self.receipt, build)
        with patch.object(audit.subprocess, "run") as run, self.assertRaisesRegex(ValueError, "clean source"):
            audit.load_build(self.receipt, "audit")
        run.assert_not_called()

    def test_incomplete_and_wrong_revision_logs_fail_closed(self):
        log = self.root / "incomplete.log"
        row = next(observations(self.revision))
        log.write_text(audit.PREFIX + json.dumps(row), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Missing"):
            audit.collect_observations(log, self.revision)
        with self.assertRaisesRegex(ValueError, "Unbound"):
            audit.collect_observations(log, "b" * 40)

    def test_fixture_identity_cannot_change_between_interactions(self):
        rows = list(observations(self.revision))
        rows[1]["fixture_sha256"] = "e" * 64
        log = self.root / "changed.log"
        log.write_text("\n".join(audit.PREFIX + json.dumps(row) for row in rows), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Fixture content changed"):
            audit.collect_observations(log, self.revision)

    def test_child_failure_retains_exit_and_exact_log_bytes(self):
        def fail(argv, **kwargs):
            if "--gentle-identity" in argv:
                return self.fake_run(argv, **kwargs)
            kwargs["stdout"].write(b"native failure\r\n")
            return subprocess.CompletedProcess(argv, 23)
        with patch.object(audit.subprocess, "run", side_effect=fail), self.assertRaisesRegex(ValueError, "23"):
            audit.run(self.receipt, self.root / "failed", "smoke", 10)
        record = json.loads((self.root / "failed/run.json").read_bytes())
        self.assertEqual(record["exit_code"], 23)
        self.assertEqual((self.root / "failed/run.log").read_bytes(), b"native failure\r\n")

    def test_timeout_is_not_success(self):
        def timeout(argv, **kwargs):
            if "--gentle-identity" in argv:
                return self.fake_run(argv, **kwargs)
            raise subprocess.TimeoutExpired(argv, 10)
        with patch.object(audit.subprocess, "run", side_effect=timeout), self.assertRaisesRegex(ValueError, "timeout"):
            audit.run(self.receipt, self.root / "failed", "smoke", 10)
        record = json.loads((self.root / "failed/run.json").read_bytes())
        self.assertIsNone(record["exit_code"])

    def test_existing_run_is_never_overwritten(self):
        output = self.root / "existing"
        output.mkdir()
        with patch.object(audit.subprocess, "run", side_effect=self.fake_run) as run, self.assertRaises(FileExistsError):
            audit.run(self.receipt, output, "smoke", 10)
        self.assertEqual(run.call_count, 1)  # Identity probe only, no benchmark.

    def enable_boundary(self):
        self.identity["workload"] = audit.BOUNDARY_WORKLOAD
        build = json.loads(self.receipt.read_bytes())
        build["workload"] = audit.BOUNDARY_WORKLOAD
        audit.write_json(self.receipt, build)

    def write_observations(self, rows, ending="\n"):
        log = self.root / "synthetic.log"
        log.write_bytes(ending.join(audit.PREFIX + json.dumps(row) for row in rows).encode())
        return log

    def test_boundary_smoke_retains_exact_case_layers_and_legacy_rows(self):
        self.enable_boundary()
        with patch.object(audit.subprocess, "run", side_effect=self.fake_run):
            result = audit.run(self.receipt, self.root / "boundary", "smoke", 10)
        self.assertTrue(result["interactive_boundary_exercised"])
        self.assertEqual(result["workload"], audit.BOUNDARY_WORKLOAD)
        self.assertFalse(result["native_gui_measured"])
        self.assertEqual(len(result["observations"]), 90)
        boundary = [row for row in result["observations"] if (row["length_bp"], row["feature_count"]) == (250_000, 5_000)]
        self.assertEqual(len(boundary), 9)
        for old, new in zip(observations(self.revision), result["observations"][:81]):
            self.assertEqual(old, {key: new[key] for key in old})
        self.assertEqual(len((self.root / "boundary/work.tsv").read_text().splitlines()), 91)

    def test_boundary_cannot_be_omitted_duplicated_or_claimed_by_legacy_binary(self):
        self.enable_boundary()
        rows = list(observations(self.revision, audit.BOUNDARY_WORKLOAD))
        for invalid in (rows[:81], rows[:-1] + [rows[0]]):
            with self.subTest(count=len(invalid)), self.assertRaisesRegex(ValueError, "Missing or duplicated"):
                audit.collect_observations(self.write_observations(invalid), self.revision, audit.BOUNDARY_WORKLOAD)
        self.identity.pop("workload")
        with patch.object(audit.subprocess, "run", side_effect=self.fake_run), self.assertRaisesRegex(ValueError, "Binary workload"):
            audit.run(self.receipt, self.root / "not-launched", "smoke", 10)
        self.assertFalse((self.root / "not-launched").exists())

    def test_unknown_workload_is_rejected(self):
        self.identity["workload"] = "future_unknown"
        with patch.object(audit.subprocess, "run", side_effect=self.fake_run), self.assertRaisesRegex(ValueError, "Unsupported DNA latency workload"):
            audit.load_build(self.receipt, "smoke")

    def test_boundary_layers_are_required_and_not_invented_zeros(self):
        for mutation in ("missing", "disabled", "empty", "negative", "wrong_feature_count", "stale_workload"):
            rows = list(observations(self.revision, audit.BOUNDARY_WORKLOAD))
            row = rows[-1]
            if mutation == "missing":
                row["layers_after"] = None
            elif mutation == "disabled":
                row["layers_after"]["orfs"]["enabled"] = False
            elif mutation == "empty":
                row["layers_after"]["restriction_groups"]["viewport_eligible"] = 0
            elif mutation == "negative":
                row["layers_before"]["gc_bins"]["total"] = -1
            elif mutation == "wrong_feature_count":
                row["layers_before"]["annotation_features"] = 4_999
            else:
                row["workload"] = audit.LEGACY_WORKLOAD
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                audit.collect_observations(self.write_observations(rows), self.revision, audit.BOUNDARY_WORKLOAD)

    def test_boundary_observations_accept_lf_and_crlf_without_changing_raw_hashes(self):
        rows = list(observations(self.revision, audit.BOUNDARY_WORKLOAD))
        lf = self.write_observations(rows)
        lf_hash = audit.sha256(lf)
        self.assertEqual(audit.collect_observations(lf, self.revision, audit.BOUNDARY_WORKLOAD), rows)
        crlf = self.write_observations(rows, "\r\n")
        self.assertNotEqual(audit.sha256(crlf), lf_hash)
        self.assertEqual(audit.collect_observations(crlf, self.revision, audit.BOUNDARY_WORKLOAD), rows)

    def test_statistics_bind_every_case_not_only_the_total_count(self):
        rows = list(observations(self.revision, audit.BOUNDARY_WORKLOAD))
        cases = {
            f"dna_feature_latency/{row['length_bp']}bp_{row['feature_count']}features_{row['fixture_sha256'][:12]}/{case}"
            for row in rows for case in audit.CASES
        }
        self.assertEqual(len(cases), 130)
        for index, case in enumerate(sorted(cases)):
            path = self.root / "criterion" / str(index) / "new"
            path.mkdir(parents=True)
            audit.write_json(path / "benchmark.json", {"full_id": case})
            audit.write_json(path / "estimates.json", {"mean": {"point_estimate": 1}})
        statistics = audit.collect_statistics(self.root / "criterion", rows)
        self.assertEqual(len(statistics), 130)
        path = self.root / "criterion/0/new/benchmark.json"
        for invalid in (sorted(cases)[1], "dna_feature_latency/wrong_fixture/constructor"):
            audit.write_json(path, {"full_id": invalid})
            with self.assertRaisesRegex(ValueError, "130 Criterion estimates"):
                audit.collect_statistics(self.root / "criterion", rows)
        path.unlink()
        with self.assertRaisesRegex(ValueError, "Missing Criterion case identity"):
            audit.collect_statistics(self.root / "criterion", rows)


if __name__ == "__main__":
    unittest.main()
