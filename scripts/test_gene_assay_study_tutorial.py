"""Source checks and optional supplied-binary replay; never builds or captures screens."""

import json
import os
from pathlib import Path
import re
import tempfile
import unittest

from prepare_gene_assay_study_tutorial import ROOT, prepare


class StudyTutorialTests(unittest.TestCase):
    def test_catalog_and_checkpoint_contract(self):
        source = json.loads((ROOT / "docs/tutorial/sources/04-08_gene_assay_study_gui.json").read_text())
        self.assertEqual(source["catalog"]["status"], "manual/hybrid")
        self.assertNotIn("generated_chapter", source)
        reviews = json.loads((ROOT / "docs/tutorial/review_manifest.json").read_text())
        entry = next(row for row in reviews["entries"] if row["tutorial_id"] == source["id"])
        self.assertEqual(entry["tutorial_kind"], "guided_walkthrough")
        for path in (ROOT / "docs/tutorial/sources").glob("*.json"):
            other = json.loads(path.read_text())
            if other.get("id") != source["id"] and "catalog" in other:
                self.assertNotEqual(other["catalog"]["order"], source["catalog"]["order"])
        guide = ROOT / source["catalog"]["path"]
        text = guide.read_text()
        for checkpoint in range(1, 7):
            self.assertIn(f"Checkpoint G{checkpoint}", text)
        for target in re.findall(r"\]\(([^)]+)\)", text):
            if "://" not in target:
                self.assertTrue((guide.parent / target.split("#")[0]).exists(), target)
        for term in ("minus-strand", "not order", "Pending", "zero-based half-open",
                     "not a promised successful study", "not a tested provider"):
            self.assertIn(term, text)

    @unittest.skipUnless(os.environ.get("GENTLE_TUTORIAL_BIN_DIR"), "supply built CLI for real replay")
    def test_real_cli_preparation_preserves_pending_study_and_candidate_provenance(self):
        binary = Path(os.environ["GENTLE_TUTORIAL_BIN_DIR"]) / ("gentle_cli.exe" if os.name == "nt" else "gentle_cli")
        with tempfile.TemporaryDirectory() as temp:
            output = prepare(binary, Path(temp) / "study")
            receipt = json.loads((output / "preparation-receipt.json").read_text())
            self.assertFalse(receipt["study_executed"])
            self.assertFalse(receipt["gui_accepted"])
            self.assertTrue(all(call["exit_code"] == 0 for call in receipt["calls"]))
            panel = json.loads((output / "patz1_sybr_juc_panel.json").read_text())
            self.assertEqual(panel["transcript_count"], 3)
            self.assertEqual(panel["strand"], "-")
            self.assertEqual(panel["selected_assay_count"], 3)
            self.assertTrue(all(row["genomic_confirmation_status"] == "not_run" for row in panel["specificity_followups"]))
            canonical = json.loads((output / "cli-dossier/canonical-report.json").read_text())
            self.assertFalse(canonical["complete"])
            self.assertEqual(canonical["pending_gene_count"], 1)
            self.assertEqual(canonical["genes"][0]["handoffs"], [])
            with self.assertRaises(FileExistsError):
                prepare(binary, output)


if __name__ == "__main__":
    unittest.main()
