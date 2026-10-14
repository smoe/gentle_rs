"""Verify retained terminal CI observations, not an execution of GENtle."""

import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
BASE = Path("docs/audits/glen_followups_20261008/candidate_de9df8af")
SOURCE = "de9df8af0a56923df01a7492f61005b1fe079d8c"
HASHES = {
    "run.api-original.json": "7f282e734ae670391679e659e9d36fa5ceb436b36b8a70c3859e31758ea60071",
    "jobs.api-original.json": "03b9987d24d1d81ada9f1b705627be083feb187d6c8661d78d76b015cfb0f5ea",
    "default_stack_from_logs.json": "d21e394e1543e86734ec061cdee91a02404b68bfbbfd5e5b0bad9b741238475f",
}


def retained_hashes():
    return {BASE / name: digest for name, digest in HASHES.items()}


def record(name):
    return json.loads((ROOT / BASE / name).read_bytes())


class FrozenCiTests(unittest.TestCase):
    def test_original_bytes_and_exact_membership(self):
        for path, digest in retained_hashes().items():
            with self.subTest(path=path):
                self.assertEqual(hashlib.sha256((ROOT / path).read_bytes()).hexdigest(), digest)
        actual = {path.relative_to(ROOT) for path in (ROOT / BASE).rglob("*") if path.is_file()}
        self.assertEqual(actual, set(retained_hashes()) | {BASE / "README.md"})

    def test_run_is_terminal_failed_and_source_bound(self):
        run = record("run.api-original.json")
        self.assertEqual(run["id"], 38599805415)
        self.assertEqual(run["head_sha"], SOURCE)
        self.assertEqual(run["head_branch"], "main")
        self.assertEqual(run["event"], "push")
        self.assertEqual(run["run_attempt"], 1)
        self.assertEqual(run["status"], "completed")
        self.assertEqual(run["conclusion"], "failure")

    def test_primary_failure_and_cascade_do_not_become_passes(self):
        payload = record("jobs.api-original.json")
        self.assertEqual(payload["total_count"], 10)
        jobs = {job["name"]: job for job in payload["jobs"]}
        self.assertEqual(len(jobs), 10)
        self.assertTrue(all(job["status"] == "completed" for job in jobs.values()))
        for platform, job_id in (("Linux", 115975878354), ("Windows", 115975818934)):
            job = jobs[platform]
            self.assertEqual(job["id"], job_id)
            self.assertEqual(job["head_sha"], SOURCE)
            self.assertEqual(job["conclusion"], "failure")
            steps = {step["name"]: step for step in job["steps"]}
            failed = [step["name"] for step in job["steps"] if step["conclusion"] == "failure"]
            self.assertEqual(failed, ["Allele and reporter fixtures without ambient stack override"])
            self.assertEqual(steps["Cargo check"]["conclusion"], "success")
            self.assertEqual(steps["Full test suite"]["conclusion"], "skipped")
            self.assertEqual(steps["Render snapshot tests"]["conclusion"], "skipped")
            self.assertEqual(steps["Cargo check (all features)"]["conclusion"], "skipped")
        self.assertEqual(jobs["CI / selected platform"]["conclusion"], "failure")

    def test_headless_safety_and_optional_skips_keep_their_scope(self):
        jobs = {job["name"]: job for job in record("jobs.api-original.json")["jobs"]}
        for name in ("Release candidate safety checks", "Linux / headless CLI and MCP"):
            self.assertEqual(jobs[name]["conclusion"], "success")
        for name in ("macOS", "08.04 synthetic GUI contracts (not online or scientific acceptance)",
                     "GUI scroll CPU smoke (not native latency)",
                     "macOS CLI for local authenticated audit (no credentials in CI)"):
            self.assertEqual(jobs[name]["conclusion"], "skipped")

    def test_removed_stack_override_is_not_an_acceptance_verdict(self):
        observation = record("default_stack_from_logs.json")
        self.assertEqual(observation["source_revision"], SOURCE)
        self.assertEqual(observation["run_id"], 38599805415)
        self.assertEqual(observation["run_attempt"], 1)
        self.assertIs(observation["default_stack_accepted"], False)
        self.assertEqual(len(observation["records"]), 2)
        self.assertEqual({row["platform"] for row in observation["records"]}, {"linux", "windows"})
        for row in observation["records"]:
            self.assertEqual(row["result"], "failed_stack_overflow")
            self.assertEqual(row["selected_tests"], 11)
            runner = row["runner"]
            self.assertEqual(runner["schema"], "gentle.default_stack_test_runner.v1")
            self.assertEqual(runner["command"][1:], ["dbsnp_", "--test-threads=1"])
            self.assertEqual(runner["inherited_rust_min_stack"], "16777216")
            self.assertIsNone(runner["child_rust_min_stack"])


if __name__ == "__main__":
    unittest.main()
