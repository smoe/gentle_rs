"""Executable companion for the hand-crafted synthetic tutorial 06.07.

Recreation and provenance: docs/tutorial/inputs/README.md. The independent
four-word oracle is literal teaching data, not a second production solver.
Actual CLI tests require GENTLE_TUTORIAL_BIN_DIR; they authorize only this
synthetic fixture, never a user's DNA. No builds, network, model or GUI capture.
"""

import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
GUIDE = ROOT / "docs/tutorial/06-07_synthetic_sequence_design.md"
SOURCE = ROOT / "docs/tutorial/sources/06-07_synthetic_sequence_design.json"
REQUEST = ROOT / "docs/tutorial/inputs/synthetic_sequence_design.json"
BIN_DIR = os.environ.get("GENTLE_TUTORIAL_BIN_DIR")
OUTPUT = "ATGTTTAAGTAA"


def fenced_request(text):
    return json.loads(re.search(r"```json\r?\n(.*?)\r?\n```", text, re.S).group(1))


class SequenceDesignTutorialSourceTests(unittest.TestCase):
    def test_discovery_and_links(self):
        source = json.loads(SOURCE.read_text(encoding="utf-8"))
        self.assertEqual(source["schema"], "gentle.tutorial_source.v4")
        self.assertNotIn("generated_chapter", source)
        self.assertEqual((source["catalog"]["group"],
                          source["catalog"]["group_position"]), ("06", 7))
        catalog = json.loads((ROOT / "docs/tutorial/catalog.json").read_text(encoding="utf-8"))
        entry, = [row for row in catalog["entries"] if row["id"] == source["id"]]
        self.assertEqual(entry["title"], source["title"])
        self.assertEqual(entry["notes"], source["catalog"]["notes"])
        self.assertEqual(entry["network"], "offline")
        text = GUIDE.read_text(encoding="utf-8")
        self.assertEqual(text.splitlines()[0], f'# {source["title"]}')
        for target in re.findall(r"\]\(([^)]+)\)", text):
            if "://" not in target:
                self.assertTrue((GUIDE.parent / target.split("#")[0]).exists(), target)

    def test_documented_request_is_exact_in_lf_and_crlf(self):
        text = GUIDE.read_bytes().decode("utf-8")
        fixture = json.loads(REQUEST.read_bytes())
        for newline in ("\n", "\r\n"):
            with self.subTest(newline=repr(newline)):
                self.assertEqual(fenced_request(text.replace("\n", newline)), fixture)
                self.assertEqual(json.loads(REQUEST.read_bytes().replace(b"\n", newline.encode())),
                                 fixture)
        self.assertEqual(fixture["purpose"], "synthetic_coding_insert")
        self.assertEqual(fixture["target"]["sequence"], "ATGTTTAAATAA")
        self.assertEqual(fixture["cds"], {"start_0based": 0, "end_0based_exclusive": 12})
        self.assertEqual(fixture["protein_sequence"], "MFK")

    def test_literal_four_word_oracle_and_teaching_table(self):
        request = json.loads(REQUEST.read_bytes())
        self.assertEqual(request["avoid_motifs"], [{"pattern": "TTTAAA", "strand": "both"}])
        self.assertEqual("TTTAAA".translate(str.maketrans("ACGT", "TGCA"))[::-1], "TTTAAA")
        bounds = request["gc_content"]
        local = request["gc_window"]
        self.assertEqual((bounds["min_basis_points"], bounds["max_basis_points"]), (1666, 1667))
        self.assertEqual(local, dict(bounds, window_bp=6))
        minimum = (bounds["min_basis_points"] * 12 + 9999) // 10000
        maximum = bounds["max_basis_points"] * 12 // 10000
        self.assertEqual((minimum, maximum), (2, 2))
        self.assertEqual(((local["min_basis_points"] * 6 + 9999) // 10000,
                          local["max_basis_points"] * 6 // 10000), (1, 1))
        # All legal words in this tiny frozen-start/stop space, independent of GENtle.
        words = ["ATGTTTAAATAA", "ATGTTCAAATAA", OUTPUT, "ATGTTCAAGTAA"]
        residues = {"ATG": "M", "TTT": "F", "TTC": "F", "AAA": "K", "AAG": "K"}
        observations = []
        for word in words:
            self.assertEqual("".join(residues[word[pos:pos + 3]] for pos in (0, 3, 6)), "MFK")
            self.assertEqual((word[:3], word[-3:]), ("ATG", "TAA"))
            gc = sum(base in "GC" for base in word)
            windows = [sum(base in "GC" for base in word[start:start + 6])
                       for start in range(7)]
            row = (gc, "TTTAAA" not in word, gc == 2, all(count == 1 for count in windows))
            observations.append(row)
            cells = [str(gc), *("yes" if value else "no" for value in row[1:])]
            self.assertIn(f'| `{word}` | ' + " | ".join(cells) + " |",
                          GUIDE.read_text(encoding="utf-8"))
        self.assertEqual(observations, [(1, False, False, False), (2, True, True, False),
                                       (2, True, True, True), (3, True, False, False)])

    def test_agent_drafts_and_review_boundaries(self):
        source = json.loads(SOURCE.read_bytes())
        text = GUIDE.read_text(encoding="utf-8")
        self.assertEqual(len(source["agent_parity"]["cases"]), 2)
        for case in source["agent_parity"]["cases"]:
            self.assertEqual(case["execution"], "ask")
            self.assertIn(case["command"], text)
        for required in ("ACTUAL_REVIEWED_DIGEST", "session", "not a digest to submit",
                         "search_claims_verified=false", "unverified_portable_preview",
                         "Do not pipe a digest automatically", "no applicable",
                         "not recommended", "outside", "Glen"):
            self.assertIn(required, text)
        self.assertIn("synthetic_sequence_design.json",
                      (REQUEST.parent / "README.md").read_text(encoding="utf-8"))


@unittest.skipUnless(BIN_DIR, "Set GENTLE_TUTORIAL_BIN_DIR for actual sequence-design CLI replay")
class SequenceDesignTutorialCliTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="gentle design's ")
        self.addCleanup(self.tmp.cleanup)
        self.cwd = Path(self.tmp.name)
        self.state = self.cwd / "tutorial.gentle.json"
        suffix = ".exe" if os.name == "nt" else ""
        self.cli = Path(BIN_DIR).resolve() / f"gentle_cli{suffix}"
        self.assertTrue(self.cli.is_file(), f"Missing explicitly configured CLI: {self.cli}")

    def invoke(self, *args, success=True):
        result = subprocess.run([str(self.cli), "--state", str(self.state), *args],
                                cwd=self.cwd, capture_output=True, timeout=90)
        if not success:
            self.assertNotEqual(result.returncode, 0, result.stdout.decode(errors="replace"))
            return result
        self.assertEqual(result.returncode, 0, result.stderr.decode(errors="replace"))
        return json.loads(result.stdout)

    def plan(self, request, name):
        path = self.cwd / f"{name} request.json"
        path.write_text(json.dumps(request), encoding="utf-8")
        report_path = self.cwd / f"{name} preview.json"
        result = self.invoke("sequence-design", "plan", f"@{path}", "--path", str(report_path))
        report = json.loads(report_path.read_bytes())
        self.assertEqual(result["result"]["dna_sequence_design"], report)
        return path, report_path, report

    def test_both_strategies_direct_and_shared_shell_preview_then_reviewed_apply(self):
        for strategy in ("conflict_directed", "full_enumeration"):
            with self.subTest(strategy=strategy):
                self.state = self.cwd / f"{strategy}.gentle.json"
                # Inline sentinel DNA is synthetic, unrelated to the redesign request.
                self.invoke("op", json.dumps({"CreateSequenceFromText": {
                    "sequence_text": "ACGT", "output_id": "keep_me", "name": None,
                    "circular": False}}))
                before = self.state.read_bytes()
                initial = json.loads(before)
                request = json.loads(REQUEST.read_bytes())
                request["search_strategy"] = strategy
                path, report_path, report = self.plan(request, strategy)
                self.assertEqual(self.state.read_bytes(), before, "Preview changed project bytes")
                self.assertEqual(report["status"], "feasible")
                self.assertEqual(report["output_sequence"], OUTPUT)
                self.assertEqual(report["edits"], [{"position_0based": 8, "before": "A", "after": "G"}])
                self.assertEqual(report["initial_matches"], [{"motif_index": 0,
                    "interval": {"start_0based": 3, "end_0based_exclusive": 9}, "strand": "both"}])
                self.assertEqual(report["gc_content"]["output"]["gc_bases"], 2)
                self.assertEqual([row["gc_bases"] for row in report["gc_window"]["input"]],
                                 [1, 1, 1, 0, 0, 0, 0])
                rows = report["gc_window"]["output"]
                self.assertEqual(len(rows), 7)
                for start, row in enumerate(rows):
                    self.assertEqual(row, {"interval": {"start_0based": start,
                        "end_0based_exclusive": start + 6}, "gc_bases": 1, "satisfies_bounds": True})
                self.assertTrue(report["optimization_complete"])
                self.assertTrue(report["minimum_edits_proven"])
                self.assertTrue(report["approval_digest"].startswith("sha256:"))
                shell_result = self.invoke("shell", "sequence-design plan " +
                                           json.dumps("@" + path.as_posix()))
                self.assertEqual(shell_result["result"]["dna_sequence_design"], report)
                self.assertEqual(self.state.read_bytes(), before)
                for args in ((), ("--approve", "sha256:not-reviewed")):
                    self.invoke("sequence-design", "apply", f"@{report_path}", *args, success=False)
                    self.assertEqual(self.state.read_bytes(), before)
                applied = self.invoke("sequence-design", "apply", f"@{report_path}",
                                      "--approve", report["approval_digest"])
                receipt = applied["result"]["dna_sequence_design_receipt"]
                self.assertTrue(receipt["output_constraints_verified"])
                self.assertFalse(receipt["search_claims_verified"])
                state = json.loads(self.state.read_bytes())
                self.assertEqual(state["sequences"]["keep_me"], initial["sequences"]["keep_me"])
                product = state["sequences"][request["output_seq_id"]]
                self.assertEqual(bytes(product["seq"]["seq"]).decode("ascii"), OUTPUT)
                record = state["metadata"]["dna_sequence_design:" + request["output_seq_id"]]
                self.assertEqual(record["submitted_proposal"], {
                    "verification": "unverified_portable_preview", "report": report})
                self.assertEqual(record["receipt"], receipt)
                saved = self.state.read_bytes()
                self.invoke("sequence-design", "apply", f"@{report_path}",
                            "--approve", report["approval_digest"], success=False)
                self.assertEqual(self.state.read_bytes(), saved, "Duplicate apply changed state")
                history = self.invoke("shell", "history status")
                self.assertEqual(history["undo_count"], 0)
                self.invoke("shell", "history undo", success=False)
                self.assertEqual(self.state.read_bytes(), saved,
                                 "A fresh CLI process must not invent undo history")

    def test_exhausted_search_has_no_approval_or_created_output(self):
        request = json.loads(REQUEST.read_bytes())
        request["max_evaluations"] = 1
        for strategy in ("conflict_directed", "full_enumeration"):
            with self.subTest(strategy=strategy):
                request["search_strategy"] = strategy
                _, path, report = self.plan(request, "exhausted " + strategy)
                self.assertEqual(report["status"], "search_exhausted")
                self.assertIsNone(report["output_sequence"])
                self.assertIsNone(report["approval_digest"])
                self.assertIsNone(report["gc_window"]["output"])
                self.invoke("sequence-design", "apply", f"@{path}",
                            "--approve", "sha256:not-reviewed", success=False)
                self.assertFalse(self.state.exists())


if __name__ == "__main__":
    unittest.main()
