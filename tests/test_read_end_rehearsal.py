"""Synthetic input/source checks and opt-in actual CLI replay, never a local build."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET

from scripts import prepare_read_end_rehearsal as rehearsal

ROOT = Path(__file__).resolve().parents[1]
BIN = os.environ.get("GENTLE_READ_END_BIN_DIR")


class ReadEndRehearsalInputs(unittest.TestCase):
    def test_oriented_synthetic_originals_are_independent_of_sam_reverse_complement(self):
        for strand in ("+", "-"):
            with self.subTest(strand=strand), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary).resolve() / "original read's α"
                rehearsal.prepare(root, oriented=True, strand=strand)
                originals = json.loads((root / "original-sequences.json").read_bytes())
                sam = [line.split("\t") for line in (root / "synthetic.sam").read_text().splitlines() if not line.startswith("@")]
                self.assertEqual({row[1] for row in sam}, {"0", "16"})
                for row in sam:
                    expected = originals[row[0]]
                    if row[1] == "16":
                        expected = expected.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                    self.assertEqual(row[9], expected)
                    self.assertTrue(originals[row[0]].endswith("T" * 18))
                model = json.loads((root / "models.json").read_bytes())
                self.assertTrue(all(m["strand"] == strand for m in model["transcripts"]))
                request = json.loads((root / "sam-request.json").read_bytes())
                manifest = json.loads((root / "units.json").read_bytes())
                self.assertEqual(manifest["sam_sha256"], rehearsal.bound(root / "synthetic.sam")["sha256"])
                (root / "saved-sam/artifacts").mkdir(parents=True)
                rehearsal.write(root / "saved-sam/report.json", {"status": "completed", "request_sha256": "synthetic", "content_sha256": "synthetic"})
                rehearsal.write(root / "saved-sam/artifacts/analysis.json", {"analysis": {
                    "context_sha256": "sha256:" + "c" * 64,
                    "units": [{"unit": {"source_sha256": request["sam"]["sha256"], "unit_id": u["unit_id"]},
                               "unclipped_axis_length": u["unclipped_axis_length"], "rows": u["rows"]} for u in manifest["units"]]}})
                rehearsal.bind_analysis(root)
                rehearsal.bind_originals(root)
                command = json.loads((root / "oriented-build-command.json").read_bytes())
                bases = command["request"]["samples"][0]["original_bases"]
                self.assertEqual(bases, rehearsal.bound(root / "original-bases.json"))
                document = json.loads(Path(bases["path"]).read_bytes())
                self.assertEqual(document["orientation"], "original_sequencing_5_to_3")
                self.assertEqual({row["read_id"]: row["sequence"] for row in document["reads"]}, originals)
                raw = Path(bases["path"]).read_bytes()
                self.assertEqual(json.loads(raw), json.loads(raw.replace(b"\n", b"\r\n")))
                self.assertNotEqual(hashlib.sha256(raw).hexdigest(), hashlib.sha256(raw.replace(b"\n", b"\r\n")).hexdigest())
                with self.assertRaises(FileExistsError):
                    rehearsal.bind_originals(root)

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

    def test_native_view_uses_shared_renderer_and_service_without_file_io(self):
        source = (ROOT / "src/app/rna_panel_ui/read_ends.rs").read_text()
        production = source.split("#[cfg(test)]", 1)[0]
        for required in ("show_rows", "render_lane", "submit_parsed", "Command::Inspect", "Command::Export", "cancel_requested", "selected_transcript_compatible_units", "gene_assigned_units"):
            self.assertIn(required, production)
        for forbidden in ("std::fs", "File::open", "project_read_ends", "geometry::", "Command::new", "std::thread"):
            self.assertNotIn(forbidden, production)
        self.assertIn("state.active_sets_report", production)
        self.assertIn("page.report_sha256", production)

    def test_tutorial_teaches_shared_synthetic_route_and_keeps_old_evidence_historical(self):
        text = (ROOT / "docs/tutorial/07-01_rna_read_batch_gene_support_cli.md").read_text()
        for required in ("Start With Saved Read Ends", "prepare_read_end_rehearsal.py", "rna-reads read-ends", "two denominators", "right censored", "no_bound_basecalls", "Historical Python Oracle", "historical captures do not certify"):
            self.assertIn(required, text)
        source = json.loads((ROOT / "docs/tutorial/sources/07-01_rna_read_batch_gene_support_cli.json").read_bytes())
        catalog = json.loads((ROOT / "docs/tutorial/catalog.json").read_bytes())
        entry, = [row for row in catalog["entries"] if row["id"] == source["id"]]
        self.assertEqual(entry["notes"], source["catalog"]["notes"])


@unittest.skipUnless(BIN, "Set GENTLE_READ_END_BIN_DIR for source-matched real CLI acceptance")
class ReadEndRehearsalProcess(unittest.TestCase):
    def test_actual_oriented_cli_both_transcript_strands_literal_runs_and_exports(self):
        binary = Path(BIN) / ("gentle_cli.exe" if os.name == "nt" else "gentle_cli")
        self.assertTrue(binary.is_file())
        for strand in ("+", "-"):
            with self.subTest(strand=strand), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary).resolve() / "oriented read's α"
                rehearsal.prepare(root, oriented=True, strand=strand)
                state = root / "not-created-project.json"
                def run(*args, success=True):
                    result = subprocess.run([str(binary), "--state", str(state), *args], capture_output=True, timeout=60)
                    self.assertEqual(result.returncode == 0, success, result.stderr.decode(errors="replace"))
                    return result
                run("rna-reads", "transcript-structures-from-sam", "@" + str(root / "sam-request.json"))
                rehearsal.bind_analysis(root)
                rehearsal.bind_originals(root)
                run("rna-reads", "read-ends", "@" + str(root / "oriented-build-command.json"))
                rehearsal.bind_pileup(root)
                run("rna-reads", "read-ends", "@" + str(root / "inspect-command.json"))
                # Inspect the published scientific page, independently of CLI response envelopes.
                report = json.loads((root / "read-ends/artifacts/pileup.json").read_bytes())
                self.assertEqual(report["schema"], "gentle.rna_read_end_pileup.v2")
                self.assertEqual(report["orientation_summary"]["same_direction_units"], 1)
                self.assertEqual(report["orientation_summary"]["opposite_direction_units"], 1)
                self.assertEqual(report["orientation_summary"]["units_with_bound_bases"], 2)
                for lane in report["lanes"]:
                    for h in lane["hypotheses"]:
                        runs = h["oriented"]["variants"][0]["runs"]
                        self.assertEqual([(r["base"], r["length"]) for r in runs], [("A", 12), ("T", 18)])
                        if lane["unit"]["unit"]["unit_id"] == "synthetic_clipped":
                            self.assertTrue(any(not r["aligned_to_reference"] for r in runs))
                run("rna-reads", "read-ends", "@" + str(root / "export-command.json"))
                artifacts = root / "exports/artifacts"
                before = {p.name: p.read_bytes() for p in artifacts.iterdir()}
                self.assertEqual({"page.json", "lanes.tsv", "pileup.svg", "pileup.png", "pileup.pdf"}, set(before))
                svg = ET.parse(artifacts / "pileup.svg").getroot()
                labels = list(svg.iter("{http://www.w3.org/2000/svg}text"))
                text = "".join("".join(label.itertext()) for label in labels)
                width, height = float(svg.attrib["width"]), float(svg.attrib["height"])
                for label in labels:
                    self.assertLessEqual(float(label.attrib["x"]) + float(label.attrib["textLength"]), width - 10)
                    self.assertLessEqual(float(label.attrib["y"]), height - 8)
                self.assertIn("12A", text)
                self.assertIn("18T", text)
                self.assertIn("Same direction", text)
                self.assertIn("Opposite", text)
                again = json.loads((root / "export-command.json").read_bytes())
                again["output_dir"] = str(root / "exports-again")
                rehearsal.write(root / "export-again-command.json", again)
                run("rna-reads", "read-ends", "@" + str(root / "export-again-command.json"))
                self.assertEqual(before, {p.name: p.read_bytes() for p in (root / "exports-again/artifacts").iterdir()},
                                 "same-host same-source export artifacts must be deterministic")
                run("rna-reads", "read-ends", "@" + str(root / "export-command.json"), success=False)
                self.assertEqual(before, {p.name: p.read_bytes() for p in artifacts.iterdir()})
                self.assertFalse(state.exists())

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
            svg = ET.parse(root / "exports/artifacts/pileup.svg").getroot()
            width, height = float(svg.attrib["width"]), float(svg.attrib["height"])
            for text in svg.iter("{http://www.w3.org/2000/svg}text"):
                self.assertLessEqual(float(text.attrib["x"]) + float(text.attrib["textLength"]), width - 10)
                self.assertLessEqual(float(text.attrib["y"]), height - 8)
            before = {p.name: p.read_bytes() for p in (root / "exports/artifacts").iterdir()}
            run("rna-reads", "read-ends", "@" + str(root / "export-command.json"), success=False)
            self.assertEqual(before, {p.name: p.read_bytes() for p in (root / "exports/artifacts").iterdir()})
            self.assertFalse(state.exists(), "artifact-only lifecycle must not create a project")
