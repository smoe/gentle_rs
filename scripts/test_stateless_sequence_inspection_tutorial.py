import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class StatelessSequenceInspectionTutorialTests(unittest.TestCase):
    def test_fixture_span_and_manual_evidence_stay_bound(self) -> None:
        fasta = ROOT / "docs/tutorial/inputs/inline_sequence_inspection_demo.fa"
        sequence = "".join(
            line.strip() for line in fasta.read_text().splitlines() if not line.startswith(">")
        )
        self.assertEqual(len(sequence), 46)

        tutorial = (ROOT / "docs/tutorial/02-02_stateless_sequence_inspection_gui_cli.md").read_text()
        self.assertIn("`0..46`", tutorial)
        self.assertNotIn("`0..47`", tutorial)
        self.assertIn("partial / blocked", tutorial)

        evidence_path = ROOT / "docs/screenshots/stateless_sequence_inspection_gui/evidence.json"
        evidence = json.loads(evidence_path.read_text())
        self.assertEqual(evidence["input"]["length_nt"], 46)
        self.assertEqual(evidence["input"]["span_0based_half_open"], [0, 46])
        self.assertEqual(evidence["score_track_gui_blocker"]["independent_fresh_process_reproductions"], 2)
        self.assertFalse(evidence["automated_gui_acceptance"])

        review_manifest = json.loads(
            (ROOT / "docs/tutorial/review_manifest.json").read_text()
        )
        review = next(
            entry
            for entry in review_manifest["entries"]
            if entry["tutorial_id"] == "stateless_sequence_inspection_gui_cli"
        )
        self.assertEqual(review["codex_reviewed_at"], "2026-10-01")

        for capture in evidence["captures"]:
            record = capture["raw_png"]
            path = ROOT / record["path"]
            payload = path.read_bytes()
            self.assertEqual(len(payload), record["size_bytes"])
            self.assertEqual(hashlib.sha256(payload).hexdigest(), record["sha256"])


if __name__ == "__main__":
    unittest.main()
