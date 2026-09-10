"""Offline regression companions for the hand-written MCP/Feature Editor guides.

All inputs are synthetic, documented in docs/tutorial/inputs/README.md. Tests
write only temporary state/evidence and never claim GUI acceptance. Set
GENTLE_TUTORIAL_BIN_DIR to run the real binary paths; without it those tests
skip explicitly. No network, external model, references or pytest dependency.
"""

import copy
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "docs/tutorial/inputs"
SPEC = importlib.util.spec_from_file_location("mcp_roundtrip", INPUTS / "mcp_roundtrip.py")
MCP_DEMO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MCP_DEMO)
BIN_DIR = os.environ.get("GENTLE_TUTORIAL_BIN_DIR")


class WalkthroughSourceTests(unittest.TestCase):
    def test_guides_are_registered_without_becoming_generated_chapters(self):
        catalog = json.loads((ROOT / "docs/tutorial/catalog.json").read_text())
        manifest = json.loads((ROOT / "docs/tutorial/manifest.json").read_text())
        reviews = json.loads((ROOT / "docs/tutorial/review_manifest.json").read_text())
        for tutorial_id, stem in [
            ("mcp_offline_roundtrip", "01-03_mcp_offline_roundtrip"),
            ("feature_editor_gui_cli", "02-05_feature_editor_gui_cli"),
        ]:
            source = json.loads((ROOT / f"docs/tutorial/sources/{stem}.json").read_text())
            self.assertEqual(source["id"], tutorial_id)
            self.assertEqual(source["catalog"]["source"], "hand_written_markdown")
            self.assertNotIn("executable", source)
            entries = [entry for entry in catalog["entries"] if entry["id"] == tutorial_id]
            self.assertEqual(len(entries), 1)
            self.assertNotEqual(entries[0]["type"], "executable_chapter")
            self.assertNotIn(tutorial_id, [c["id"] for c in manifest["chapters"]])
            self.assertIn(tutorial_id, [r["tutorial_id"] for r in reviews["entries"]])
            guide = ROOT / source["catalog"]["path"]
            for target in re.findall(r"\]\(([^)]+)\)", guide.read_text()):
                if "://" not in target:
                    self.assertTrue((guide.parent / target.split("#")[0]).exists(), target)

    def test_fixture_recreation_and_provenance(self):
        text = (INPUTS / "feature_editor_demo.gb").read_text()
        bases = re.sub(r"[^acgt]", "", text.split("ORIGIN\n")[1].split("//")[0])
        self.assertEqual(bases, "acgt" * 30)
        self.assertIn("All features are invented", text)
        for name in ["feature_editor_demo.gb", "mcp_restriction_operation.json", "mcp_roundtrip.py"]:
            self.assertIn(name, (INPUTS / "README.md").read_text())

    def test_comparison_excludes_only_execution_identity(self):
        report = {"report_id": "a", "op_id": "one", "run_id": "x",
                  "generated_at_unix_ms": 1, "rows": [{"recognition_start_0based": 0}],
                  "scan_topology": "linear", "target_label": "demo"}
        other = dict(report, report_id="b", op_id="two", run_id="y", generated_at_unix_ms=2)
        self.assertEqual(MCP_DEMO.comparable_report(report), MCP_DEMO.comparable_report(other))
        changed = copy.deepcopy(other)
        changed["rows"][0]["recognition_start_0based"] = 1
        self.assertNotEqual(MCP_DEMO.comparable_report(report), MCP_DEMO.comparable_report(changed))
        self.assertEqual(report["op_id"], "one")

    def test_mcp_framing_and_truncated_response(self):
        messages = [{"jsonrpc": "2.0", "id": 1, "result": {"label": "\u00e4"}},
                    {"jsonrpc": "2.0", "id": 2, "result": {}}]
        raw = MCP_DEMO.encode_frames(messages)
        self.assertEqual(MCP_DEMO.decode_frames(raw), messages)
        for invalid in [raw[:-1], b"status text\n", b"Content-Length: -1\r\n\r\n{}"]:
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    MCP_DEMO.decode_frames(invalid)

    def test_tool_errors_are_not_successful_process_verdicts(self):
        for message in [{"error": {"code": -32602}},
                        {"result": {"isError": True, "content": [{"text": "refused"}]}}]:
            with self.assertRaises(RuntimeError):
                MCP_DEMO.tool_payload(message)


@unittest.skipUnless(BIN_DIR, "Set GENTLE_TUTORIAL_BIN_DIR for real MCP/CLI tutorial replay")
class McpWalkthroughTests(unittest.TestCase):
    def test_discovery_confirmation_and_complete_report_parity(self):
        suffix = ".exe" if os.name == "nt" else ""
        cli, mcp = [Path(BIN_DIR).resolve() / f"gentle_{name}{suffix}" for name in ["cli", "mcp"]]
        with tempfile.TemporaryDirectory(prefix="gentle-mcp-tutorial-") as tmp:
            run_dir = Path(tmp) / "run"
            report = MCP_DEMO.run_demo(cli, mcp, run_dir)
            self.assertEqual(report["matched_site_count"], 2)
            self.assertEqual(report["source_sequence_length_bp"], 12)
            self.assertTrue(all(row["forward_strand"] for row in report["rows"]))
            self.assertEqual(json.loads((run_dir / "summary.json").read_text())["status"], "pass")
            before = (run_dir / "summary.json").read_bytes()
            with self.assertRaises(FileExistsError):
                MCP_DEMO.run_demo(cli, mcp, run_dir)
            self.assertEqual((run_dir / "summary.json").read_bytes(), before)


@unittest.skipUnless(BIN_DIR, "Set GENTLE_TUTORIAL_BIN_DIR for real Feature Editor command replay")
class FeatureEditorWalkthroughTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="gentle-editor-tutorial-")
        self.addCleanup(self.tmp.cleanup)
        self.state = Path(self.tmp.name) / "editor.gentle.json"
        suffix = ".exe" if os.name == "nt" else ""
        self.cli = Path(BIN_DIR).resolve() / f"gentle_cli{suffix}"
        self.invoke("op", json.dumps({"LoadFile": {
            "path": str(INPUTS / "feature_editor_demo.gb"), "as_id": "editor_demo"}}))
        self.original = self.sequence_record()

    def invoke(self, *args, expect_success=True):
        result = subprocess.run([self.cli, "--state", self.state, *args], cwd=ROOT,
                                capture_output=True, timeout=60)
        if expect_success:
            self.assertEqual(result.returncode, 0, result.stderr.decode(errors="replace"))
            return json.loads(result.stdout)
        self.assertNotEqual(result.returncode, 0, result.stdout.decode(errors="replace"))
        return result

    def shell(self, command, **kwargs):
        return self.invoke("shell", command, **kwargs)

    def sequence_record(self):
        return json.loads(self.state.read_text())["sequences"]["editor_demo"]["seq"]

    def preview(self, command):
        before = self.state.read_bytes()
        report = self.shell(command + " --dry-run")["report"]
        self.assertEqual(self.state.read_bytes(), before, "Preview mutated the saved project")
        self.assertTrue(report["dry_run"])
        self.assertFalse(report["applied"])
        return report

    def assert_rows(self, intervals):
        result = self.shell("features query editor_demo --sort feature_id --include-qualifiers")
        self.assertEqual(result["sequence_length_bp"], 120)
        self.assertEqual(result["total_feature_count"], len(intervals))
        self.assertEqual([r["feature_id"] for r in result["rows"]], list(range(len(intervals))))
        self.assertEqual([(r["start_0based"], r["end_0based_exclusive"]) for r in result["rows"]], intervals)
        return result["rows"]

    def test_preview_location_create_split_merge_delete_and_reopen(self):
        self.assert_rows([(10, 100), (20, 40), (60, 80)])
        command = "features edit-location editor_demo 1 --start-1based 21 --end-1based-inclusive 45"
        preview = self.preview(command)
        lock = preview["before_feature_fingerprint_sha256"]
        self.assertEqual(preview["after"]["start_0based"], 20)
        self.assertEqual(preview["after"]["end_0based_exclusive"], 45)
        applied = self.shell(command + f" --expected-feature-fingerprint-sha256 {lock}")["report"]
        self.assertTrue(applied["applied"])
        unchanged = self.state.read_bytes()
        self.shell(command + f" --expected-feature-fingerprint-sha256 {lock}", expect_success=False)
        self.assertEqual(self.state.read_bytes(), unchanged, "Stale apply mutated state")

        create = ("features create editor_demo --kind misc_feature --start-1based 25 "
                  "--end-1based-inclusive 35 --strand forward "
                  "--qualifier label=review_patch --qualifier gene=DEMO")
        preview = self.preview(create)
        self.assertIsNone(preview["outcome"].get("created_feature_index"))
        candidates = preview["review_candidates"]
        self.assertEqual({c["feature_index"] for c in candidates}, {0, 1})
        self.assertTrue(all(not c["modified"] for c in candidates))
        for candidate in candidates:
            self.assertEqual({e["evidence_kind"] for e in candidate["evidence"]},
                             {"overlapping_location", "shared_identifier"})
        annotation_lock = preview["before_annotation_state_fingerprint_sha256"]
        created = self.shell(create + f" --expected-annotation-state-fingerprint-sha256 {annotation_lock}")["report"]
        self.assertEqual(created["outcome"]["created_feature_index"], 3)
        self.assert_rows([(10, 100), (20, 45), (60, 80), (24, 35)])

        split = "features split editor_demo 1 --split-before-1based 31"
        preview = self.preview(split)
        self.assertEqual(preview["outcome"]["split_at_0based"], 30)
        feature_lock = preview["outcome"]["original_feature"]["feature_fingerprint_sha256"]
        annotation_lock = preview["before_annotation_state_fingerprint_sha256"]
        self.shell(split + f" --expected-feature-fingerprint-sha256 {feature_lock}"
                   f" --expected-annotation-state-fingerprint-sha256 {annotation_lock}")
        self.assert_rows([(10, 100), (20, 30), (30, 45), (60, 80), (24, 35)])
        features = self.sequence_record()["features"]
        self.assertEqual(features[1]["qualifiers"], features[2]["qualifiers"])

        merge = "features merge editor_demo 1 2"
        preview = self.preview(merge)
        left, right = [f["feature_fingerprint_sha256"] for f in preview["outcome"]["source_features"]]
        annotation_lock = preview["before_annotation_state_fingerprint_sha256"]
        self.shell(merge + f" --expected-first-feature-fingerprint-sha256 {left}"
                   f" --expected-second-feature-fingerprint-sha256 {right}"
                   f" --expected-annotation-state-fingerprint-sha256 {annotation_lock}")
        self.assert_rows([(10, 100), (20, 45), (60, 80), (24, 35)])

        delete = "features delete editor_demo 3"
        preview = self.preview(delete)
        feature_lock = preview["outcome"]["deleted_feature"]["feature_fingerprint_sha256"]
        annotation_lock = preview["before_annotation_state_fingerprint_sha256"]
        self.shell(delete + f" --expected-feature-fingerprint-sha256 {feature_lock}"
                   f" --expected-annotation-state-fingerprint-sha256 {annotation_lock}")
        self.assert_rows([(10, 100), (20, 45), (60, 80)])
        final = self.sequence_record()
        self.assertEqual({k: v for k, v in final.items() if k != "features"},
                         {k: v for k, v in self.original.items() if k != "features"})
        self.assertEqual(final["features"][0], self.original["features"][0])
        self.assertEqual(final["features"][2], self.original["features"][2])
        self.assertEqual(final["features"][1]["qualifiers"], self.original["features"][1]["qualifiers"])
        reopened = Path(self.tmp.name) / "reopened.gentle.json"
        result = subprocess.run([self.cli, "--state", self.state, "save-project", reopened],
                                capture_output=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr.decode(errors="replace"))
        self.state = reopened
        self.assert_rows([(10, 100), (20, 45), (60, 80)])
        self.assertEqual(self.sequence_record(), final)

    def test_invalid_geometry_and_reverse_orientation(self):
        before = self.state.read_bytes()
        for command in [
            "features edit-location editor_demo 1 --start-1based 0 --end-1based-inclusive 45 --dry-run",
            "features split editor_demo 1 --split-before-1based 21 --dry-run",
            "features merge editor_demo 1 2 --dry-run",
        ]:
            with self.subTest(command=command):
                self.shell(command, expect_success=False)
                self.assertEqual(self.state.read_bytes(), before)
        preview = self.preview("features edit-location editor_demo 2 --start-1based 61 --end-1based-inclusive 82")
        self.assertEqual(preview["after"]["strand"], "reverse")
        self.assertEqual(preview["after"]["five_prime_position_1based"], 82)
        self.assertEqual(preview["after"]["three_prime_position_1based"], 61)


if __name__ == "__main__":
    unittest.main()
