"""Source checks and optional supplied-binary replay; never builds or captures screens."""

import json
import os
from pathlib import Path
import re
import struct
import tempfile
import unittest

from prepare_real_patz1_tutorial import ROOT, FIXTURE, prepare, sha


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
            self.assertIn(term, " ".join(text.split()).replace("**", ""))

    def test_authentic_patz1_inputs_are_pinned_not_synthetic(self):
        manifest = json.loads((FIXTURE / "manifest.json").read_text())
        for name, expected in manifest["files"].items():
            self.assertEqual(sha(FIXTURE / name), expected)
        entry = json.loads((FIXTURE / "ensembl_entry.json").read_text())
        self.assertEqual(entry["gene_id"], "ENSG00000100105")
        self.assertEqual(len(entry["transcripts"]), 13)
        self.assertEqual(entry["strand"], -1)
        self.assertEqual(entry["sequence_length"], 20802)

    def test_public_gui_checkpoint_evidence_is_hash_bound(self):
        root = ROOT / "docs/screenshots/gene_assay_study_gui"
        evidence = json.loads((root / "evidence.json").read_text())
        self.assertEqual(evidence["schema"], "gentle.manual_gui_tutorial_evidence.v1")
        self.assertEqual(evidence["status"], "manual_hybrid_pass")
        self.assertFalse(evidence["automated_gui_acceptance"])
        self.assertFalse(evidence["inner_agent_invoked"])
        self.assertEqual(
            evidence["producer"]["revision"],
            "3c1c32bcceacbe7d327c09f83ba57dcca2682b06",
        )
        self.assertEqual(evidence["result"]["completion_status"], "partial")
        self.assertEqual(evidence["result"]["transcript_records"], 13)
        self.assertEqual(evidence["result"]["exact_mature_cdna_classes"], 13)
        self.assertEqual(evidence["result"]["selected_assays"], 7)
        self.assertEqual(evidence["result"]["unresolved_class_pairs"], 9)
        self.assertEqual(evidence["result"]["genomic_specificity_assessments"], 0)
        self.assertEqual(evidence["result"]["whole_transcriptome_specificity"], "not_run")

        self.assertEqual([row["id"] for row in evidence["captures"]], [f"G{i}" for i in range(1, 7)])
        for row in evidence["captures"]:
            for field in ("raw_png", "semantic_snapshot", "context_svg"):
                path = root / row[field]
                self.assertTrue(path.is_file(), path)
                self.assertEqual(sha(path), row[f"{field}_sha256"], path)
            raw = (root / row["raw_png"]).read_bytes()
            self.assertEqual(raw[:8], b"\x89PNG\r\n\x1a\n")
            self.assertEqual(struct.unpack(">II", raw[16:24]), (1920, 1080))
            snapshot = json.loads((root / row["semantic_snapshot"]).read_text())
            self.assertEqual(snapshot["schema"], "gentle.gui_semantic_snapshot.v2")
            self.assertGreater(len(snapshot["items"]), 0)
            context = (root / row["context_svg"]).read_text()
            self.assertIn(row["raw_png"], context)

    @unittest.skipUnless(os.environ.get("GENTLE_TUTORIAL_BIN_DIR"), "supply built CLI for real replay")
    def test_real_cli_preparation_preserves_pending_study_and_candidate_provenance(self):
        binary = Path(os.environ["GENTLE_TUTORIAL_BIN_DIR"]) / ("gentle_cli.exe" if os.name == "nt" else "gentle_cli")
        with tempfile.TemporaryDirectory() as temp:
            output = prepare(binary, Path(temp) / "study")
            receipt = json.loads((output / "preparation-receipt.json").read_text())
            self.assertFalse(receipt["study_executed"])
            self.assertFalse(receipt["gui_accepted"])
            self.assertTrue(all(call["exit_code"] == 0 for call in receipt["calls"]))
            self.assertFalse(receipt["synthetic"])
            self.assertFalse(receipt["primer_design_executed"])
            locus = json.loads((output / "locus.report.json").read_text())
            self.assertEqual(locus["isoform_evidence"]["splicing"]["transcript_count"], 13)
            self.assertEqual(len(locus["transcript_presentation"]["records"]), 17)
            operation = json.loads((output / "primer-discrimination.operation.json").read_text())["DesignTranscriptAssayPanel"]
            self.assertEqual(operation["objective"], "minimal_discrimination_panel")
            self.assertNotIn("junction_evidence_paths", operation)
            publication = json.loads((output / "publication.request.json").read_text())
            self.assertEqual(publication["genes"][0]["status"], "pending")
            self.assertFalse((output / "primer-panel.json").exists())
            with self.assertRaises(FileExistsError):
                prepare(binary, output)


if __name__ == "__main__":
    unittest.main()
