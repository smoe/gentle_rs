"""Offline maintenance regressions using synthetic, disposable repositories.

All documents below are hand-crafted in memory, written only to temporary
directories, and recreated by this suite. No project, biological fixture,
generated artifact, Cargo invocation or network access is needed.
"""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from scripts import maintenance_chore as chore


ROOT = Path(__file__).resolve().parents[1]


class RoadmapBudgetTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="GENtle O'Brien \u00dc roadmap ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "docs").mkdir()
        self.roadmap = self.root / "docs/roadmap.md"

    def write_roadmap(self, text):
        self.roadmap.write_bytes(text.encode("utf-8"))

    def test_word_limit_is_inclusive_and_overflow_is_only_a_warning(self):
        for count in (0, 1, 1999, 2000, 2001):
            with self.subTest(words=count):
                self.write_roadmap(" ".join(["word"] * count))
                result = chore.check_roadmap_word_count(self.root)
                self.assertEqual(result.status, "warn" if count > 2000 else "ok")
                self.assertIn(f"has {count} whitespace-separated words", result.evidence[0])
                self.assertEqual(bool(result.findings), count > 2000)
                self.assertEqual(bool(result.minimal_fix), count > 2000)

    def test_wrapping_tabs_and_lf_crlf_do_not_change_the_word_budget(self):
        for newline in ("\n", "\r\n"):
            for separator in (" ", "\t", newline):
                with self.subTest(newline=newline, separator=separator):
                    self.write_roadmap(separator.join(["word"] * 2000) + newline)
                    result = chore.check_roadmap_word_count(self.root)
                    self.assertEqual(result.status, "ok")
                    self.assertIn("has 2000 whitespace-separated words", result.evidence[0])

    def test_compact_but_dense_roadmap_warns_independently_of_line_count(self):
        self.write_roadmap(" ".join(["word"] * 2001) + "\n" * 200)
        self.assertEqual(chore.check_roadmap_line_count(self.root).status, "ok")
        self.assertEqual(chore.check_roadmap_word_count(self.root).status, "warn")

    def test_line_limits_remain_unchanged(self):
        for count in (199, 200, 400, 401):
            with self.subTest(lines=count):
                self.write_roadmap("word\n" * count)
                result = chore.check_roadmap_line_count(self.root)
                self.assertEqual(result.status, "ok" if 200 <= count <= 400 else "fail")

    def prepare_synthetic_repo(self, newline, lines):
        documents = {
            "docs/architecture.md": "docs/decisions.md docs/CHANGELOG.md\n",
            "docs/CHANGELOG.md": "# Synthetic changelog\n",
            "docs/decisions.md": "# Synthetic decisions\n",
            "docs/maintenance_chore_plan.md": "# Synthetic maintenance plan\n",
            "AGENTS.md": (
                "docs/maintenance_chore_plan.md docs/CHANGELOG.md\n"
                "Do not add a `Done` entry\n"
            ),
        }
        for name, text in documents.items():
            (self.root / name).write_bytes(text.replace("\n", newline).encode("utf-8"))
        text = "CHANGELOG.md decisions.md Do not add a `Done` entry Move completed work directly "
        text += " ".join(["word"] * 2001) + newline * lines
        self.write_roadmap(text)
        subprocess.run(["git", "init", "--quiet", str(self.root)], check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def run_session_close(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts/maintenance_chore.py"),
             "--repo-root", str(self.root), "session-close", "--format", "json",
             "--plan", "synthetic roadmap budget regression"],
            cwd=self.root, capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.stderr, "")
        return result, json.loads(result.stdout)

    def test_real_session_close_warns_without_failing_for_lf_and_crlf(self):
        for newline in ("\n", "\r\n"):
            with self.subTest(newline=newline):
                self.prepare_synthetic_repo(newline, 200)
                result, report = self.run_session_close()
                self.assertEqual(result.returncode, 0, report)
                checks = {item["name"]: item for item in report["checks"]}
                self.assertEqual(checks["roadmap word count"]["status"], "warn")
                self.assertEqual(checks["roadmap line count"]["status"], "ok")
                self.assertEqual(report["summary"]["fail"], 0)
                self.assertIn("do not drop open requirements",
                              checks["roadmap word count"]["minimal_fix"][0])

    def test_real_session_close_still_fails_other_violations(self):
        self.prepare_synthetic_repo("\n", 199)
        result, report = self.run_session_close()
        self.assertEqual(result.returncode, 1)
        checks = {item["name"]: item for item in report["checks"]}
        self.assertEqual(checks["roadmap word count"]["status"], "warn")
        self.assertEqual(checks["roadmap line count"]["status"], "fail")


if __name__ == "__main__":
    unittest.main()
