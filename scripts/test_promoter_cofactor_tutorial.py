"""Synthetic tutorial metadata checks and opt-in shared-CLI/real-Parquet replay.

No Cargo builds, network, screenshots or publication. Runtime replay requires
both GENTLE_TUTORIAL_BIN_DIR and GENTLE_TEST_DUCKDB; configured bad paths fail.
Fixture origin/recreation: test_files/fixtures/promoter_cofactors/README.md.
"""

import hashlib
from datetime import date
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from promoter_cofactor_tutorial import ROOT, FIXTURE, prepare, replay

SOURCE = ROOT / "docs/tutorial/sources/08-14_promoter_cofactor_browser.json"
BIN_DIR = os.environ.get("GENTLE_TUTORIAL_BIN_DIR")
DUCKDB = os.environ.get("GENTLE_TEST_DUCKDB")


def validate_human_review(review):
    """Pending and genuine dated sign-offs are both valid; partial claims are not."""
    reviewed_at = review.get("human_reviewed_at")
    reviewer = review.get("human_reviewer")
    if reviewed_at is None and reviewer is None:
        return
    if not isinstance(reviewer, str) or not reviewer.strip():
        raise ValueError("A dated human review requires a named reviewer")
    if not isinstance(reviewed_at, str) or date.fromisoformat(reviewed_at).isoformat() != reviewed_at:
        raise ValueError("Human review date must be YYYY-MM-DD")


class PromoterCofactorTutorial(unittest.TestCase):
    def test_registered_source_and_teaching_links(self):
        source = json.loads(SOURCE.read_text(encoding="utf-8"))
        catalog = source["catalog"]
        self.assertEqual(source["schema"], "gentle.tutorial_source.v4")
        self.assertEqual((catalog["group"], catalog["group_position"]), ("08", 14))
        self.assertIn("agent_users", catalog["audiences"])
        self.assertIn("ui open tutorial-guide promoter_cofactor_browser", catalog["notes"])
        generated = json.loads((ROOT / "docs/tutorial/catalog.json").read_text(encoding="utf-8"))
        entry = next(e for e in generated["entries"] if e["id"] == source["id"])
        self.assertEqual(entry["path"], catalog["path"])
        self.assertEqual(entry["notes"], catalog["notes"])
        guide = ROOT / catalog["path"]
        text = guide.read_text(encoding="utf-8")
        self.assertEqual(text.splitlines()[0], "# " + source["title"])
        for word in ("17 bp", "1/2", "source_record", "4.30329", "MA0653.1",
                     "not a binding-probability ratio", "unavailable", "half-open"):
            self.assertIn(word, text)
        for relative in re.findall(r"\]\(([^)]+)\)", text):
            self.assertTrue((guide.parent / relative.split("#")[0]).is_file(), relative)
        reviews = json.loads((ROOT / "docs/tutorial/review_manifest.json").read_text(encoding="utf-8"))
        review = next(r for r in reviews["entries"] if r["tutorial_id"] == source["id"])
        validate_human_review(review)

    def test_human_review_allows_pending_or_complete_signoff(self):
        validate_human_review({"human_reviewed_at": None, "human_reviewer": None})
        validate_human_review({"human_reviewed_at": "2026-09-16", "human_reviewer": "synthetic-test-reviewer"})
        for date_value, reviewer in ((None, "reviewer"), ("2026-09-16", None),
                                     ("not-a-date", "reviewer"), ("2026-09-16", " ")):
            with self.subTest(date=date_value, reviewer=reviewer), self.assertRaises(ValueError):
                validate_human_review({"human_reviewed_at": date_value, "human_reviewer": reviewer})

    def test_fixture_template_declares_scope_and_score_family(self):
        template = json.loads((FIXTURE / "manifest.template.json").read_text(encoding="utf-8"))
        self.assertFalse(template["complete_genome_scan"])
        self.assertEqual(len(template["distance_bands"]), 6)
        self.assertEqual(template["source_score_floor"], -1)
        self.assertEqual(template["positive_threshold"], 0)
        self.assertEqual(template["score_configuration"]["score_mode"], "log2_relative_risk")
        self.assertEqual(template["files"], [])
        self.assertIn("manifest.template.json", (ROOT / "src/promoter_cofactors_tests.rs").read_text(encoding="utf-8"))

    def test_existing_package_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            sentinel = path / "scientific-output.txt"
            sentinel.write_bytes(b"keep original")
            with self.assertRaises(FileExistsError):
                prepare(path, Path("never-invoke-this"))
            self.assertEqual(sentinel.read_bytes(), b"keep original")

    @unittest.skipUnless(BIN_DIR and DUCKDB, "Set GENTLE_TUTORIAL_BIN_DIR and GENTLE_TEST_DUCKDB for real CLI/Parquet")
    def test_shared_cli_replay_captures_exact_source_report(self):
        cli = (Path(BIN_DIR) / ("gentle_cli.exe" if os.name == "nt" else "gentle_cli")).resolve(strict=True)
        duckdb = Path(DUCKDB).resolve(strict=True)
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "path with spaces and 'quotes' and \u00f6"
            with patch("promoter_cofactor_tutorial.subprocess.run", wraps=subprocess.run) as run:
                saved = replay(output, cli, duckdb)
                for call in run.call_args_list:
                    self.assertEqual(call.kwargs.get("encoding"), "utf-8")
            evidence = saved["evidence"][0]
            canonical = json.dumps(evidence["source_record"], ensure_ascii=False,
                                   separators=(",", ":"), sort_keys=True).encode()
            self.assertEqual(evidence["source_sha256"], "sha256:" + hashlib.sha256(canonical).hexdigest())
            self.assertEqual(saved["interval"]["strand"], "plus")
            self.assertIsNone(saved.get("local_projection"))
            self.assertNotIn("max_signal_value", evidence)
            self.assertEqual(evidence["associated_gene_ids"], ["GENE-A", "GENE-B"])
            receipt = json.loads((output / "replay.json").read_text(encoding="utf-8"))
            self.assertEqual(receipt["feature_sequences"], ["demo_plus"])
            for seq_id, strand in (("demo_plus", "plus"),):
                self.assertTrue((output / f"{seq_id}.preview.state.json").is_file())
                report = json.loads((output / f"{seq_id}.applied.json").read_text(encoding="utf-8"))
                self.assertEqual(report["feature_materialization"]["projection"]["local_strand"], strand)
                self.assertTrue(report["feature_materialization"]["curation"]["applied"])
            for command in receipt["commands"]:
                self.assertEqual(command[1:3], ["--state", str(output / "tutorial.state.json")])


if __name__ == "__main__":
    unittest.main()
