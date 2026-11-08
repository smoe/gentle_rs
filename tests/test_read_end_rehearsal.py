"""Synthetic input/source checks and opt-in actual CLI replay, never a local build."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from scripts import prepare_read_end_rehearsal as rehearsal

ROOT = Path(__file__).resolve().parents[1]
BIN = os.environ.get("GENTLE_READ_END_BIN_DIR")


class ReadEndRehearsalInputs(unittest.TestCase):
    def test_bound_synthetic_models_annotation_and_manifest_under_both_newlines(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve() / "rehearsal's α"
            rehearsal.prepare(root)
            request = json.loads((root / "sam-request.json").read_bytes())
            self.assertEqual(request["required_flank_bases"], 2)
            for name in ("sam", "manifest", "dictionary"):
                file = request[name]
                raw = Path(file["path"]).read_bytes()
                self.assertEqual(file["sha256"], "sha256:" + hashlib.sha256(raw).hexdigest())
                self.assertEqual(file["bytes"], len(raw))
            model = json.loads((root / "models.json").read_bytes())
            annotation = json.loads((root / "annotation.json").read_bytes())
            self.assertEqual(model["schema"], "gentle.rna_read_end_transcript_models.v2")
            for ending in ("\n", "\r\n"):
                self.assertEqual(json.loads(json.dumps(model, indent=2).replace("\n", ending)), model)
            self.assertEqual(model["annotation_sha256"], rehearsal.bound(root / "annotation.json")["sha256"])
            for a, m in zip(annotation["transcripts"], model["transcripts"]):
                self.assertEqual(a["genomic_exons"], m["genomic_exons_0based_half_open"])
            self.assertEqual(model["transcripts"][0]["genomic_stop_codon_0based_half_open"], [[220, 223]])
            manifest = json.loads((root / "units.json").read_bytes())
            self.assertTrue(all(not u["search"]["views_complete"] and not u["search"]["search_complete"] for u in manifest["units"]))
            with self.assertRaises(FileExistsError):
                rehearsal.prepare(root)

    def test_command_drafts_keep_exact_bindings_and_all_view_inputs_explicit(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve() / "fixture"
            rehearsal.prepare(root)
            (root / "saved-sam").mkdir()
            rehearsal.write(root / "saved-sam/report.json", {"status": "completed", "request_sha256": "synthetic request", "content_sha256": "synthetic content"})
            rehearsal.bind_analysis(root)
            command = json.loads((root / "build-command.json").read_bytes())
            self.assertEqual(command["request"]["transcript_id"], "SYN_TX_A")
            self.assertEqual(command["request"]["samples"][0]["composition"]["orientation"], {"status": "undeclared"})
            self.assertIsNone(command["request"]["samples"][0]["tail"])
            (root / "read-ends").mkdir()
            rehearsal.write(root / "read-ends/report.json", {"synthetic": True})
            rehearsal.bind_pileup(root)
            inspect = json.loads((root / "inspect-command.json").read_bytes())
            self.assertEqual(inspect["view"], rehearsal.view())
            self.assertEqual(inspect["report"], rehearsal.bound(root / "read-ends/report.json"))
            export = json.loads((root / "export-command.json").read_bytes())
            self.assertEqual(export["formats"], ["json", "tsv", "svg", "png", "pdf"])

    def test_shared_route_and_policy_are_named_without_frontend_geometry(self):
        policy = (ROOT / "src/engine/operation_lifecycle.rs").read_text()
        self.assertIn("Operation::RnaReadEnds", policy)
        self.assertIn("command.requires_confirmation()", policy)
        managed = (ROOT / "src/command_execution/rna_artifacts.rs").read_text()
        for required in ("Policy::ReadEnds", "Command::Cancel", "completion_binding"):
            self.assertIn(required, managed)
        mcp = (ROOT / "src/mcp_server.rs").read_text()
        self.assertIn('"rna-reads read-ends"', mcp)
        service = (ROOT / "src/rna_read_ends.rs").read_text()
        for required in ("inspect_observed_transcript_structures_from_sam", "read_sam_bound_json", "publish_read_end_bundle", "core::project_read_ends", "pileup::build"):
            self.assertIn(required, service)


@unittest.skipUnless(BIN, "Set GENTLE_READ_END_BIN_DIR for source-matched real CLI acceptance")
class ReadEndRehearsalProcess(unittest.TestCase):
    def test_actual_cli_prepare_build_inspect_export_and_no_clobber(self):
        binary = Path(BIN) / ("gentle_cli.exe" if os.name == "nt" else "gentle_cli")
        self.assertTrue(binary.is_file())
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve() / "rehearsal's α"
            rehearsal.prepare(root)
            state = root / "not-created-project.json"
            def run(*args, success=True):
                result = subprocess.run([str(binary), "--state", str(state), *args], capture_output=True, timeout=60)
                self.assertEqual(result.returncode == 0, success, result.stderr.decode(errors="replace"))
                return result
            run("rna-reads", "transcript-structures-from-sam", "@" + str(root / "sam-request.json"))
            rehearsal.bind_analysis(root)
            run("rna-reads", "read-ends", "@" + str(root / "build-command.json"))
            rehearsal.bind_pileup(root)
            run("rna-reads", "read-ends", "@" + str(root / "inspect-command.json"))
            run("rna-reads", "read-ends", "@" + str(root / "export-command.json"))
            before = {p.name: p.read_bytes() for p in (root / "exports/artifacts").iterdir()}
            run("rna-reads", "read-ends", "@" + str(root / "export-command.json"), success=False)
            self.assertEqual(before, {p.name: p.read_bytes() for p in (root / "exports/artifacts").iterdir()})
            self.assertFalse(state.exists(), "artifact-only lifecycle must not create a project")
