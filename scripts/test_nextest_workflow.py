"""Offline checks for local nextest setup; never build or execute Rust tests."""

from pathlib import Path
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[1]


class NextestWorkflowTests(unittest.TestCase):
    def test_recording_is_explicit_and_retention_is_bounded(self):
        config = tomllib.loads((ROOT / ".config/nextest-user.toml").read_text())
        self.assertIs(config["experimental"]["record"], True)
        self.assertIs(config["record"]["enabled"], True)
        self.assertEqual(config["record"]["max-records"], 100)
        self.assertEqual(config["record"]["max-total-size"], "1GB")
        self.assertEqual(config["record"]["max-age"], "30d")

    def test_repository_profile_does_not_hide_failures_or_oversubscribe(self):
        config = tomllib.loads((ROOT / ".config/nextest.toml").read_text())
        profile = config["profile"]["default"]
        self.assertEqual(profile["test-threads"], 1)
        self.assertIs(profile["fail-fast"], False)
        self.assertEqual(profile["retries"], 0)
        self.assertNotIn("overrides", profile)
        self.assertNotIn("default-filter", profile)
        self.assertNotIn("experimental", config)
        self.assertNotIn("record", config)

    def test_reruns_are_bracketed_by_full_runs_and_doctests(self):
        docs = (ROOT / "docs/testing.md").read_text()
        section = docs.split("### Local Failure Collection And Recorded Reruns", 1)[1].split(
            "## 2. Engine tests", 1
        )[0]
        self.assertIn("cargo test --locked -j1 --no-fail-fast", section)
        command = "cargo nextest run --locked --build-jobs 1 --no-fail-fast"
        self.assertEqual(section.count(command + "\n"), 2)
        self.assertIn(command + " --rerun latest", section)
        self.assertIn("cargo test --locked -j1 --doc", section)
        self.assertIn("--user-config-file .config/nextest-user.toml", section)
        self.assertIn("one process per test", section)
        self.assertIn("not a\nfull-suite verdict", section)


if __name__ == "__main__":
    unittest.main()
